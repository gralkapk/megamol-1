#pragma once

#include <chrono>
#include <fstream>
#include <string>

#include "Screenshots.h"

#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libswscale/swscale.h>

#include <glm/glm.hpp>

namespace megamol::frontend {
struct StreamContext {
    StreamContext() = default;
    ~StreamContext() {
        if (enc_ctx != nullptr) {
            avcodec_free_context(&enc_ctx);
        }
        if (dec_ctx != nullptr) {
            avcodec_free_context(&dec_ctx);
        }
        if (rgbpic != nullptr) {
            av_frame_free(&rgbpic);
        }
        if (yuvpic != nullptr) {
            av_frame_free(&yuvpic);
        }
        if (sws_ctx != nullptr) {
            sws_freeContext(sws_ctx);
        }
        if (packet != nullptr) {
            av_packet_free(&packet);
        }
        /*if (fmt_ctx != nullptr) {
            avformat_free_context(fmt_ctx);
        }*/
    }

    AVFormatContext* fmt_ctx = nullptr;
    AVFormatContext* in_fmt_ctx = nullptr;

    AVCodecContext* enc_ctx = nullptr;
    AVCodecContext* dec_ctx = nullptr;

    SwsContext* sws_ctx = nullptr;

    AVFrame* rgbpic = nullptr;
    AVFrame* yuvpic = nullptr;

    AVPacket* packet = nullptr;

    glm::ivec2 dim;
};

bool setup_video(std::string const& out_filename, glm::ivec2 const& dim, std::vector<StreamContext>& stream_ctx) {
    int ret = 0;
    AVFormatContext* ovid_fmtctx = nullptr;
    ret = avformat_alloc_output_context2(&ovid_fmtctx, nullptr, nullptr, out_filename.c_str());
    if (ret < 0)
        return false;

    stream_ctx.clear();
    stream_ctx.resize(2);

    // establish two streams: (0) video, (1) subtitle
    auto vid_stream = avformat_new_stream(ovid_fmtctx, nullptr);
    if (vid_stream == nullptr)
        return false;
    auto vid_enc = avcodec_find_encoder(AV_CODEC_ID_H264);
    if (vid_enc == nullptr)
        return false;
    auto vid_enc_ctx = avcodec_alloc_context3(vid_enc);
    if (vid_enc_ctx == nullptr)
        return false;
    vid_enc_ctx->width = dim.x;
    vid_enc_ctx->height = dim.y;
    vid_enc_ctx->pix_fmt = AV_PIX_FMT_YUV420P;
    vid_enc_ctx->time_base = AVRational{1, 30};
    if (ovid_fmtctx->oformat->flags & AVFMT_GLOBALHEADER)
        vid_enc_ctx->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;
    ret = avcodec_open2(vid_enc_ctx, vid_enc, nullptr);
    if (ret < 0)
        return false;
    ret = avcodec_parameters_from_context(vid_stream->codecpar, vid_enc_ctx);
    if (ret < 0)
        return false;
    vid_stream->time_base = vid_enc_ctx->time_base;
    stream_ctx[0].enc_ctx = vid_enc_ctx;

    auto sub_header = std::string(
        "[Script Info]\r\n; Script generated by FFmpeg/Lavc58.134.100\r\nScriptType: v4.00+\r\nPlayResX: "
        "384\r\nPlayResY: 288\r\nScaledBorderAndShadow: yes\r\n\r\n[V4+ Styles]\r\nFormat: Name, Fontname, "
        "Fontsize, PrimaryColour, SecondaryColour, OutlineColour, BackColour, Bold, Italic, Underline, "
        "StrikeOut, ScaleX, ScaleY, Spacing, Angle, BorderStyle, Outline, Shadow, Alignment, MarginL, MarginR, "
        "MarginV, Encoding\r\nStyle: "
        "Default,Arial,16,&Hffffff,&Hffffff,&H0,&H0,0,0,0,0,100,100,0,0,1,1,0,2,10,10,10,0\r\n\r\n[Events]"
        "\r\nFormat: Layer, Start, End, Style, Name, MarginL, MarginR, MarginV, Effect, Text");

    auto sub_stream = avformat_new_stream(ovid_fmtctx, nullptr);
    if (sub_stream == nullptr)
        return false;
    auto sub_enc = avcodec_find_encoder(AV_CODEC_ID_SUBRIP);
    if (sub_enc == nullptr)
        return false;
    auto sub_enc_ctx = avcodec_alloc_context3(sub_enc);
    if (sub_enc_ctx == nullptr)
        return false;
    sub_enc_ctx->time_base = AVRational{1, 30};
    sub_enc_ctx->subtitle_header = (uint8_t*)sub_header.c_str();
    sub_enc_ctx->subtitle_header_size = strlen(sub_header.c_str());
    ret = avcodec_open2(sub_enc_ctx, sub_enc, nullptr);
    if (ret < 0)
        return false;
    sub_stream->time_base = sub_enc_ctx->time_base;
    stream_ctx[1].enc_ctx = sub_enc_ctx;

    av_dump_format(ovid_fmtctx, 0, out_filename.c_str(), 1);

    avio_open(&ovid_fmtctx->pb, out_filename.c_str(), AVIO_FLAG_WRITE);
    avformat_write_header(ovid_fmtctx, nullptr);

    stream_ctx[0].sws_ctx = sws_getContext(
        dim.x, dim.y, AV_PIX_FMT_RGB24, dim.x, dim.y, AV_PIX_FMT_YUV420P, SWS_FAST_BILINEAR, nullptr, nullptr, nullptr);

    stream_ctx[0].rgbpic = av_frame_alloc();
    stream_ctx[0].rgbpic->format = AV_PIX_FMT_RGB24;
    stream_ctx[0].rgbpic->width = dim.x;
    stream_ctx[0].rgbpic->height = dim.y;
    ret = av_frame_get_buffer(stream_ctx[0].rgbpic, 1);
    if (ret < 0)
        return false;

    stream_ctx[0].yuvpic = av_frame_alloc();
    stream_ctx[0].yuvpic->format = AV_PIX_FMT_YUV420P;
    stream_ctx[0].yuvpic->width = dim.x;
    stream_ctx[0].yuvpic->height = dim.y;
    ret = av_frame_get_buffer(stream_ctx[0].yuvpic, 1);
    if (ret < 0)
        return false;

    stream_ctx[0].dim = dim;
    stream_ctx[0].packet = av_packet_alloc();
    stream_ctx[0].fmt_ctx = ovid_fmtctx;
    stream_ctx[1].dim = dim;
    stream_ctx[1].packet = av_packet_alloc();
    stream_ctx[1].fmt_ctx = ovid_fmtctx;

    return false;
}

void rgb2yuv(StreamContext& config) {
    sws_scale(config.sws_ctx, config.rgbpic->data, config.rgbpic->linesize, 0, config.dim.y, config.yuvpic->data,
        config.yuvpic->linesize);
}

void flipRGB(StreamContext& config, megamol::frontend_resources::ScreenshotImageData const& image) {
    for (int y = 0; y < config.dim.y; ++y) {
        for (int x = 0; x < config.dim.x; ++x) {
            config.rgbpic->data[0][y * config.rgbpic->linesize[0] + 3 * x] = image.flipped_rows[y][x].r;
            config.rgbpic->data[0][y * config.rgbpic->linesize[0] + 3 * x + 1] = image.flipped_rows[y][x].g;
            config.rgbpic->data[0][y * config.rgbpic->linesize[0] + 3 * x + 2] = image.flipped_rows[y][x].b;
        }
    }
}

void write_srt_entry(
    std::ofstream& srt, uint64_t id, std::string const& base_ts, std::string const& next_ts, std::string const& msg) {
    srt << std::to_string(id) << std::endl;
    srt << base_ts << " --> " << next_ts << std::endl;
    srt << msg << std::endl;
}

void encodeFrame(StreamContext& config) {
    av_packet_unref(config.packet);
    auto enc_ret = avcodec_send_frame(config.enc_ctx, config.yuvpic);
    while (enc_ret >= 0) {
        enc_ret = avcodec_receive_packet(config.enc_ctx, config.packet);
        if (enc_ret == AVERROR(EAGAIN) || enc_ret == AVERROR_EOF)
            break;
        config.packet->stream_index = 0;
        enc_ret = av_interleaved_write_frame(config.fmt_ctx, config.packet);
    }
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

bool setup_subtitles(std::string const& filename, std::vector<StreamContext>& stream_ctx) {
    //AVFormatContext* isub_fmtctx = nullptr;
    int ret = 0;
    ret = avformat_open_input(&stream_ctx[1].in_fmt_ctx, filename.c_str(), nullptr, nullptr);
    if (ret < 0)
        return false;

    ret = avformat_find_stream_info(stream_ctx[1].in_fmt_ctx, nullptr);
    if (ret < 0)
        return false;

    if (stream_ctx[1].in_fmt_ctx->nb_streams > 1)
        return false;

    auto stream = stream_ctx[1].in_fmt_ctx->streams[0];
    auto dec = avcodec_find_decoder(stream->codecpar->codec_id);
    auto codec_ctx = avcodec_alloc_context3(dec);
    ret = avcodec_parameters_to_context(codec_ctx, stream->codecpar);
    if (ret < 0)
        return false;
    codec_ctx->time_base = AVRational{1, 30};
    ret = avcodec_open2(codec_ctx, dec, nullptr);
    if (ret < 0)
        return false;
    stream_ctx[1].dec_ctx = codec_ctx;

    //for (int i = 0; i < stream_ctx[1].in_fmt_ctx->nb_streams; ++i) {
    //    AVStream* stream = stream_ctx[1].in_fmt_ctx->streams[i];
    //    AVCodec const* dec = avcodec_find_decoder(stream->codecpar->codec_id);
    //    AVCodecContext* codec_ctx = avcodec_alloc_context3(dec);
    //    avcodec_parameters_to_context(codec_ctx, stream->codecpar);
    //    /*if (codec_ctx->codec_type == AVMEDIA_TYPE_VIDEO
    //            || codec_ctx->codec_type == AVMEDIA_TYPE_AUDIO) {
    //            if (codec_ctx->codec_type == AVMEDIA_TYPE_VIDEO)
    //                    codec_ctx->framerate = av_guess_frame_rate(ivid_fmtctx, stream, nullptr);
    //            FFEXEC(avcodec_open2(codec_ctx, dec, nullptr));
    //    }*/
    //    codec_ctx->time_base = AVRational{1, 30};
    //    avcodec_open2(codec_ctx, dec, nullptr);
    //    stream_sub_ctx[i].dec_ctx = codec_ctx;
    //}

    //av_dump_format(isub_fmtctx, 0, filename.c_str(), 0);
}

void encode_sub(std::vector<StreamContext>& stream_ctx) {
    int ret = 0;
    while (1) {
        av_packet_unref(stream_ctx[1].packet);
        ret = av_read_frame(stream_ctx[1].in_fmt_ctx, stream_ctx[1].packet);
        //int got_sub = 0;
        //AVSubtitle* sub = new AVSubtitle;
        //avcodec_decode_subtitle2(stream_ctx[1].dec_ctx, sub, &got_sub, stream_ctx[1].packet);
        /*if (!got_sub)
            break;*/
        stream_ctx[1].packet->stream_index = 1;
        //av_packet_rescale_ts(packet, stream_sub_ctx[0].dec_ctx->time_base, stream_ctx[1].enc_ctx->time_base);
        /*packet->time_base = AVRational{ 1,30 };
        packet->pts = sub_counter;
        sub_counter += 10;*/
        ret = av_interleaved_write_frame(stream_ctx[1].fmt_ctx, stream_ctx[1].packet);
        /*uint8_t* buf = new uint8_t[2048 * 2048];
        ret = avcodec_encode_subtitle(stream_ctx[1].enc_ctx, buf, 2048 * 2048, sub);*/
        if (ret < 0)
            break;
    }
}
} // namespace megamol::frontend
