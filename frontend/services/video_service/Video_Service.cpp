/*
 * Video_Service.cpp
 *
 * Copyright (C) 2021 by MegaMol Team
 * Alle Rechte vorbehalten.
 */

// search/replace Template_Service with your class name
// you should also delete the FAQ comments in these template files after you read and understood them
#include "Video_Service.hpp"

#include "mmcore/MegaMolGraph.h"

#include <chrono>
#include <fstream>
#include <sstream>
#include <string>

#include "Screenshots.h"

#include "glad/gl.h"


//#define VPL_ENABLED


// local logging wrapper for your convenience until central MegaMol logger established
#include "mmcore/utility/log/Log.h"

#ifdef VPL_ENABLED
#include "vpl/mfx.h"
#define VPLVERSION(major, minor) (major << 16 | minor)
#define ALIGN16(value) (((value + 15) >> 4) << 4)
#endif

extern "C" {
#include "libavcodec/avcodec.h"
#include "libavformat/avformat.h"
#include "libavutil/avutil.h"
}

using fps_30 = std::chrono::duration<double, std::ratio<1, 30>>;

std::string convert_to_timestamp(fps_30 const& time) {
    auto t_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time);

    auto t_s = std::chrono::duration_cast<std::chrono::seconds>(t_ms);
    t_ms -= t_s;

    auto t_m = std::chrono::duration_cast<std::chrono::minutes>(t_s);
    t_s -= t_m;

    auto t_h = std::chrono::duration_cast<std::chrono::hours>(t_m);
    t_m -= t_h;

    char buf[13];
    sprintf_s(buf, 13, "%02d:%02d:%02d,%03d", t_h.count(), t_m.count(), t_s.count(), t_ms.count());


    return std::string(buf);
}

void write_srt_entry(
    std::ofstream& srt, uint64_t id, std::string const& base_ts, std::string const& next_ts, std::string const& msg) {
    srt << std::to_string(id) << std::endl;
    srt << base_ts << " --> " << next_ts << std::endl;
    srt << msg << std::endl;
}

static const std::string service_name = "Video_Service: ";
static void log(std::string const& text) {
    const std::string msg = service_name + text;
    megamol::core::utility::log::Log::DefaultLog.WriteInfo(msg.c_str());
}

static void log_error(std::string const& text) {
    const std::string msg = service_name + text;
    megamol::core::utility::log::Log::DefaultLog.WriteError(msg.c_str());
}

static void log_warning(std::string const& text) {
    const std::string msg = service_name + text;
    megamol::core::utility::log::Log::DefaultLog.WriteWarn(msg.c_str());
}


void print_log(std::string const& msg, int ret) {
    char err_buf[AV_ERROR_MAX_STRING_SIZE];
    fprintf(stderr, "%s: %s\n", msg.c_str(), av_make_error_string(err_buf, AV_ERROR_MAX_STRING_SIZE, ret));
}


namespace megamol {
namespace frontend {

struct VPL_Setup {
#ifdef VPL_ENABLED
    mfxLoader loader;
    mfxSession session;
    mfxFrameSurface1* encSurface;
    mfxBitstream bitstream;
    std::ofstream file;
    std::ofstream srt;
    uint64_t frameID;
#endif
    AVFormatContext* oc;
    AVOutputFormat* o_fmt;
    AVCodec* o_codec;
    AVStream* ov_stream;
    AVCodecContext* enc_ctx;
    AVFrame* enc_frame;
    AVPacket* pkt;
    int64_t next_pts = 0;
};

void fill_yuv_image(AVFrame* pict, int frame_index, int width, int height) {
    int x, y, i;

    i = frame_index;

    /* Y */
    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
            pict->data[0][y * pict->linesize[0] + x] = x + y + i * 3;

    /* Cb and Cr */
    for (y = 0; y < height / 2; y++) {
        for (x = 0; x < width / 2; x++) {
            pict->data[1][y * pict->linesize[1] + x] = 128 + y + i * 2;
            pict->data[2][y * pict->linesize[2] + x] = 64 + x + i * 5;
        }
    }
}

void fill_yuv_image(AVFrame* pict, int frame_index, int width, int height,
    megamol::frontend_resources::ScreenshotImageData_I420 const& input) {
    int x, y, i;

    i = frame_index;

    /* Y */
    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
            pict->data[0][y * pict->linesize[0] + x] = input.Y[y * pict->linesize[0] + x];

    /* Cb and Cr */
    for (y = 0; y < height / 2; y++) {
        for (x = 0; x < width / 2; x++) {
            pict->data[1][y * pict->linesize[1] + x] = input.U[y * pict->linesize[1] + x];
            pict->data[2][y * pict->linesize[2] + x] = input.V[y * pict->linesize[2] + x];
        }
    }
}

AVFrame* get_video_frame(std::unique_ptr<VPL_Setup>& setup_, megamol::frontend_resources::ScreenshotImageData_I420 const& input) {
    av_frame_make_writable(setup_->enc_frame);
    //fill_yuv_image(setup_->enc_frame, setup_->next_pts, setup_->enc_ctx->width, setup_->enc_ctx->height);
    fill_yuv_image(setup_->enc_frame, setup_->next_pts, setup_->enc_ctx->width, setup_->enc_ctx->height, input);
    setup_->enc_frame->pts = setup_->next_pts++;
    return setup_->enc_frame;
}

AVFrame* get_video_frame(std::unique_ptr<VPL_Setup>& setup_) {
    av_frame_make_writable(setup_->enc_frame);
    fill_yuv_image(setup_->enc_frame, setup_->next_pts, setup_->enc_ctx->width, setup_->enc_ctx->height);
    setup_->enc_frame->pts = setup_->next_pts++;
    return setup_->enc_frame;
}

Video_Service::Video_Service() {
    // init members to default states
}

Video_Service::~Video_Service() {
    // clean up raw pointers you allocated with new, which is bad practice and nobody does
}

bool Video_Service::init(void* configPtr) {
    if (configPtr == nullptr)
        return false;

    return init(*static_cast<Config*>(configPtr));
}

bool Video_Service::init(const Config& config) {
    // initialize your service and its provided resources using config parameters
    // for now, you dont need to worry about your service beeing initialized or closed multiple times
    // init() and close() only get called once in the lifetime of each service object
    // but maybe more instances of your service will get created? this may be relevant for central resources you manage (like libraries, network connections).

    m_requestedResourcesNames = {"MegaMolGraph"};

#ifdef VPL_ENABLED
    setup_ = std::make_unique<VPL_Setup>();

    setup_->loader = MFXLoad();

    mfxConfig cfg[3];
    mfxVariant cfgVal[3];
    mfxStatus sts;

    cfg[0] = MFXCreateConfig(setup_->loader);

    cfgVal[0].Type = MFX_VARIANT_TYPE_U32;
    cfgVal[0].Data.U32 = MFX_IMPL_TYPE_SOFTWARE;

    sts = MFXSetConfigFilterProperty(cfg[0], (mfxU8*)"mfxImplDescription.Impl", cfgVal[0]);

    cfg[1] = MFXCreateConfig(setup_->loader);

    cfgVal[1].Type = MFX_VARIANT_TYPE_U32;
    cfgVal[1].Data.U32 = MFX_CODEC_HEVC;

    sts = MFXSetConfigFilterProperty(
        cfg[1], (mfxU8*)"mfxImplDescription.mfxEncoderDescription.encoder.CodecID", cfgVal[1]);

    cfg[2] = MFXCreateConfig(setup_->loader);

    cfgVal[2].Type = MFX_VARIANT_TYPE_U32;
    cfgVal[2].Data.U32 = VPLVERSION(2, 2);

    sts = MFXSetConfigFilterProperty(cfg[2], (mfxU8*)"mfxImplDescription.ApiVersion.Version", cfgVal[2]);

    sts = MFXCreateSession(setup_->loader, 0, &setup_->session);

    mfxVideoParam encodeParams = {};

    encodeParams.mfx.CodecId = MFX_CODEC_HEVC;
    encodeParams.mfx.TargetUsage = MFX_TARGETUSAGE_BALANCED;
    encodeParams.mfx.TargetKbps = 4000;
    encodeParams.mfx.RateControlMethod = MFX_RATECONTROL_VBR;
    encodeParams.mfx.FrameInfo.FrameRateExtN = 30;
    encodeParams.mfx.FrameInfo.FrameRateExtD = 1;
    encodeParams.mfx.FrameInfo.FourCC = MFX_FOURCC_I420;
    encodeParams.mfx.FrameInfo.ChromaFormat = MFX_CHROMAFORMAT_YUV420;
    encodeParams.mfx.FrameInfo.CropW = 1920;
    encodeParams.mfx.FrameInfo.CropH = 1080;
    encodeParams.mfx.FrameInfo.Width = ALIGN16(1920);
    encodeParams.mfx.FrameInfo.Height = ALIGN16(1080);

    encodeParams.IOPattern = MFX_IOPATTERN_IN_SYSTEM_MEMORY;

    sts = MFXVideoENCODE_Init(setup_->session, &encodeParams);

    sts = MFXMemory_GetSurfaceForEncode(setup_->session, &setup_->encSurface);

    setup_->bitstream = {};
    setup_->bitstream.MaxLength = 9000000;
    setup_->bitstream.Data = new mfxU8[setup_->bitstream.MaxLength]();

    setup_->file = std::ofstream("test.h265", std::ios_base::binary);
    setup_->srt = std::ofstream("test.srt");

    setup_->frameID = 0;
#endif

#ifdef DEBUG
    av_log_set_level(AV_LOG_DEBUG);
#endif //  DEBUG

    std::string out_filename = "test.mp4";

    av_register_all();
    avcodec_register_all();

    setup_ = std::make_unique<VPL_Setup>();

    auto ret = avformat_alloc_output_context2(&setup_->oc, nullptr, nullptr, out_filename.c_str());
    if (ret < 0)
        print_log("Failed to alloc context", ret);

    setup_->o_fmt = setup_->oc->oformat;

    setup_->o_codec = avcodec_find_encoder(setup_->o_fmt->video_codec);

    setup_->ov_stream = avformat_new_stream(setup_->oc, nullptr);

    setup_->ov_stream->id = setup_->oc->nb_streams - 1;

    setup_->enc_ctx = avcodec_alloc_context3(setup_->o_codec);

    setup_->enc_ctx->codec_id = setup_->o_fmt->video_codec;
    setup_->enc_ctx->bit_rate = 400000;
    setup_->enc_ctx->width = 1920;
    setup_->enc_ctx->height = 1080;
    setup_->ov_stream->time_base = AVRational{1, 30};
    setup_->enc_ctx->time_base = setup_->ov_stream->time_base;
    setup_->enc_ctx->gop_size = 12;
    setup_->enc_ctx->pix_fmt = AV_PIX_FMT_YUV420P;

    if (setup_->o_fmt->flags & AVFMT_GLOBALHEADER) {
        setup_->enc_ctx->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;
    }

    ret = avcodec_open2(setup_->enc_ctx, setup_->o_codec, nullptr);
    if (ret < 0)
        print_log("Failed to open codec", ret);

    auto srt_codec = avcodec_find_encoder(AVCodecID::AV_CODEC_ID_SUBRIP);
    auto srt_stream = avformat_new_stream(setup_->oc, nullptr);
    srt_stream->id = setup_->oc->nb_streams - 1;
    auto srt_ctx = avcodec_alloc_context3(srt_codec);
    //srt_ctx->sub_charenc = "utf-8";
    /*srt_stream->time_base = setup_->ov_stream->time_base;
    srt_ctx->time_base = srt_stream->time_base;*/
    //srt_ctx->time_base = AVRational{1, AV_TIME_BASE};
    srt_ctx->width = 1920;
    srt_ctx->height = 1080;
    srt_ctx->time_base = AVRational{1, 30};
    srt_ctx->framerate = AVRational{30, 1};
    ret = avcodec_open2(srt_ctx, srt_codec, nullptr);
    if (ret < 0)
        print_log("Failed to open codec", ret);

    setup_->enc_frame = av_frame_alloc();
    setup_->enc_frame->format = setup_->enc_ctx->pix_fmt;
    setup_->enc_frame->width = setup_->enc_ctx->width;
    setup_->enc_frame->height = setup_->enc_ctx->height;
    ret = av_frame_get_buffer(setup_->enc_frame, 0);
    if (ret < 0)
        print_log("Failed to get buffer for enc frame", ret);

    ret = avcodec_parameters_from_context(setup_->ov_stream->codecpar, setup_->enc_ctx);
    if (ret < 0)
        print_log("Failed to get copy codec par", ret);

    av_dump_format(setup_->oc, 0, out_filename.c_str(), 1);

    ret = avio_open(&setup_->oc->pb, out_filename.c_str(), AVIO_FLAG_WRITE);
    if (ret < 0)
        print_log("Failed to open output file", ret);

    ret = avformat_write_header(setup_->oc, nullptr);
    if (ret <0)
        print_log("Failed to write outout header", ret);

    setup_->pkt = av_packet_alloc();

    log("initialized successfully");
    return true;
}

void Video_Service::close() {
    // close libraries or APIs you manage
    // wrap up resources your service provides, but don not depend on outside resources to be available here
    // after this, at some point only the destructor of your service gets called

#ifdef VPL_ENABLED
    mfxSyncPoint syncp;
    auto sts = MFXVideoENCODE_EncodeFrameAsync(setup_->session, nullptr, nullptr, &setup_->bitstream, &syncp);
    if (sts == MFX_ERR_NONE && syncp) {
        // Encode output is not available on CPU until sync operation
        // completes
        sts = MFXVideoCORE_SyncOperation(setup_->session, syncp, 0);

        setup_->file.write(reinterpret_cast<char const*>(setup_->bitstream.Data + setup_->bitstream.DataOffset),
            setup_->bitstream.DataLength);
        setup_->bitstream.DataLength = 0;
    }

    if (setup_->encSurface) {
        setup_->encSurface->FrameInterface->Release(setup_->encSurface);
    }

    if (setup_->session) {
        MFXVideoENCODE_Close(setup_->session);
        MFXClose(setup_->session);
    }

    if (setup_->loader)
        MFXUnload(setup_->loader);

    delete[] setup_->bitstream.Data;

    setup_->file.close();
    setup_->srt.close();
#endif

    auto ret = av_write_trailer(setup_->oc);
    if (ret < 0)
        print_log("Failed to write trailer", ret);
}

std::vector<FrontendResource>& Video_Service::getProvidedResources() {

    return m_providedResourceReferences;
}

const std::vector<std::string> Video_Service::getRequestedResourceNames() const {
    // since this function should not change the state of the service
    // you should assign your requested resource names in init()


    return m_requestedResourcesNames;
}

void Video_Service::setRequestedResources(std::vector<FrontendResource> resources) {
    megamolgraph_ptr =
        const_cast<megamol::core::MegaMolGraph*>(&resources[0].getResource<megamol::core::MegaMolGraph>());
}

void Video_Service::updateProvidedResources() {}

void Video_Service::digestChangedRequestedResources() {}

void Video_Service::resetProvidedResources() {
    // this gets called at the end of the main loop iteration
    // since the current resources state should have been handled in this frame already
    // you may clean up resources whose state is not needed for the next iteration
    // e.g. m_keyboardEvents.clear();
    // network_traffic_buffer.reset_to_empty();
}

void Video_Service::preGraphRender() {
    // this gets called right before the graph is told to render something
    // e.g. you can start a start frame timer here

    // rendering via MegaMol View is called after this function finishes
    // in the end this calls the equivalent of ::mmcRenderView(hView, &renderContext)
    // which leads to view.Render()
}

void Video_Service::postGraphRender() {
    // the graph finished rendering and you may more stuff here
    // e.g. end frame timer
    // update window name
    // swap buffers, glClear


#ifdef VPL_ENABLED
    std::string project = megamolgraph_ptr->Convenience().SerializeGraph();

    mfxStatus sts;

    mfxSyncPoint syncp;

    GLint viewport_dims[4] = {0};
    glGetIntegerv(GL_VIEWPORT, viewport_dims);
    GLint fbWidth = viewport_dims[2];
    GLint fbHeight = viewport_dims[3];

    megamol::frontend_resources::ScreenshotImageData result;
    result.resize(static_cast<size_t>(fbWidth), static_cast<size_t>(fbHeight));

    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, fbWidth, fbHeight, GL_RGBA, GL_UNSIGNED_BYTE, result.image.data());

    megamol::frontend_resources::ScreenshotImageData_I420 input;
    input = megamol::frontend_resources::ScreenshotImageData_I420(
        megamol::frontend_resources::ScreenshotImageData_YUV(result));

    // set data in surface
    sts = setup_->encSurface->FrameInterface->Map(setup_->encSurface, MFX_MAP_WRITE);

    mfxU16 w, h, i, pitch;
    size_t bytes_read;
    mfxU8* ptr;
    mfxFrameInfo* info = &setup_->encSurface->Info;
    mfxFrameData* data = &setup_->encSurface->Data;

    /*w = info->Width;
    h = info->Height;*/
    w = input.width;
    h = input.height;

    pitch = data->Pitch;
    ptr = data->Y;
    for (i = 0; i < h; i++) {
        std::copy(input.Y.begin() + i * pitch, input.Y.begin() + (i + 1) * pitch, ptr + i * pitch);
    }

    // read chrominance (U, V)
    pitch /= 2;
    h /= 2;
    w /= 2;
    ptr = data->U;
    for (i = 0; i < h; i++) {
        std::copy(input.U.begin() + i * pitch, input.U.begin() + (i + 1) * pitch, ptr + i * pitch);
    }

    ptr = data->V;
    for (i = 0; i < h; i++) {
        std::copy(input.V.begin() + i * pitch, input.V.begin() + (i + 1) * pitch, ptr + i * pitch);
    }

    sts = setup_->encSurface->FrameInterface->Unmap(setup_->encSurface);

    sts = MFXVideoENCODE_EncodeFrameAsync(setup_->session, nullptr, setup_->encSurface, &setup_->bitstream, &syncp);

    switch (sts) {
    case MFX_ERR_NONE:
        // MFX_ERR_NONE and syncp indicate output is available
        if (syncp) {
            // Encode output is not available on CPU until sync operation
            // completes
            sts = MFXVideoCORE_SyncOperation(setup_->session, syncp, 0);

            setup_->file.write(reinterpret_cast<char const*>(setup_->bitstream.Data + setup_->bitstream.DataOffset),
                setup_->bitstream.DataLength);
            setup_->bitstream.DataLength = 0;
        }
        break;
    case MFX_ERR_NOT_ENOUGH_BUFFER:
        // This example deliberatly uses a large output buffer with immediate
        // write to disk for simplicity. Handle when frame size exceeds
        // available buffer here
        break;
    case MFX_ERR_MORE_DATA:
        // The function requires more data to generate any output
        break;
    case MFX_ERR_DEVICE_LOST:
        // For non-CPU implementations,
        // Cleanup if device is lost
        break;
    case MFX_WRN_DEVICE_BUSY:
        // For non-CPU implementations,
        // Wait a few milliseconds then try again
        break;
    default:
        printf("unknown status %d\n", sts);
        break;
    }


    /*auto test = fps_30(setup_->frameID).count();
    auto t_ms = std::chrono::duration_cast<std::chrono::milliseconds>(fps_30(setup_->frameID)).count();
    auto t_s = std::chrono::duration_cast<std::chrono::seconds>(fps_30(setup_->frameID)).count();
    auto t_m = std::chrono::duration_cast<std::chrono::minutes>(fps_30(setup_->frameID)).count();
    auto t_h = std::chrono::duration_cast<std::chrono::hours>(fps_30(setup_->frameID)).count();
    std::chrono::high_resolution_clock::time_point tp(
        std::chrono::duration_cast<std::chrono::milliseconds>(fps_30(setup_->frameID)));*/
    /* auto t_c = std::chrono::system_clock::to_time_t(tp);*/
    //std::stringstream ss;
    //ss << std::put_time(std::localtime(&t_c), "%T");

    auto base_ts = convert_to_timestamp(fps_30(setup_->frameID));
    auto next_ts = convert_to_timestamp(fps_30(setup_->frameID + 1));

    write_srt_entry(setup_->srt, setup_->frameID + 1, base_ts, next_ts, project);


    ++setup_->frameID;
#endif

    GLint viewport_dims[4] = {0};
    glGetIntegerv(GL_VIEWPORT, viewport_dims);
    GLint fbWidth = viewport_dims[2];
    GLint fbHeight = viewport_dims[3];

    megamol::frontend_resources::ScreenshotImageData result;
    result.resize(static_cast<size_t>(fbWidth), static_cast<size_t>(fbHeight));

    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, fbWidth, fbHeight, GL_RGBA, GL_UNSIGNED_BYTE, result.image.data());

    megamol::frontend_resources::ScreenshotImageData_I420 input;
    input = megamol::frontend_resources::ScreenshotImageData_I420(
        megamol::frontend_resources::ScreenshotImageData_YUV(result));

    //auto ret = avcodec_send_frame(setup_->enc_ctx, get_video_frame(setup_));
    auto ret = avcodec_send_frame(setup_->enc_ctx, get_video_frame(setup_, input));
    if (ret < 0)
        print_log("Unable to send frame", ret);

    while (ret >= 0) {
        ret = avcodec_receive_packet(setup_->enc_ctx, setup_->pkt);
        if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF)
            break;
        else if (ret < 0) {
            print_log("Unable to receive pkt", ret);
        }

        av_packet_rescale_ts(setup_->pkt, setup_->enc_ctx->time_base, setup_->ov_stream->time_base);
        setup_->pkt->stream_index = setup_->ov_stream->index;

        ret = av_interleaved_write_frame(setup_->oc, setup_->pkt);
        if (ret < 0)
            print_log("Uanble to write frame", ret);
    }

}


} // namespace frontend
} // namespace megamol
