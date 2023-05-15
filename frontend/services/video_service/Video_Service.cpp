#include "Video_Service.hpp"

#include <sstream>

#include "video_util.hpp"

#include "glad/gl.h"


#include "LuaCallbacksCollection.h"


megamol::frontend::Video_Service::Video_Service() {}


megamol::frontend::Video_Service::~Video_Service() {}


bool megamol::frontend::Video_Service::init(void* configPtr) {

    //requestedResourcesNames_ = {"RegisterLuaCallback", "MegaMolGraph"};
    requestedResourcesNames_ = {"MegaMolGraph"};

    // for test purposes
    start_video_rec("./test_out.mkv");
    srt_file_ = std::ofstream("./test_out.srt");

    image_.resize(1920, 1080);

    return true;
}


void megamol::frontend::Video_Service::close() {
    stop_video_rec("./test_out.mkv");
}


std::vector<megamol::frontend::FrontendResource>& megamol::frontend::Video_Service::getProvidedResources() {
    return providedResources_;
}


const std::vector<std::string> megamol::frontend::Video_Service::getRequestedResourceNames() const {
    return requestedResourcesNames_;
}


void megamol::frontend::Video_Service::setRequestedResources(std::vector<FrontendResource> resources) {
    mmgraph_ptr = const_cast<megamol::core::MegaMolGraph*>(&resources[0].getResource<megamol::core::MegaMolGraph>());
}


void megamol::frontend::Video_Service::updateProvidedResources() {}


void megamol::frontend::Video_Service::digestChangedRequestedResources() {}


void megamol::frontend::Video_Service::resetProvidedResources() {}


void megamol::frontend::Video_Service::preGraphRender() {
    // handle start and stop
    // need lua callbacks for that
}


std::string parameter_diff(std::string const& lhs, std::string const& rhs) {
    std::istringstream lhs_stream(lhs);
    std::istringstream rhs_stream(rhs);
    std::stringstream diff;
    for (std::string lhs_line, rhs_line; std::getline(lhs_stream, lhs_line), std::getline(rhs_stream, rhs_line);) {
        if (lhs_line != rhs_line) {
            diff << rhs_line << "\n";
        }
    }

    return diff.str();
}


void megamol::frontend::Video_Service::postGraphRender() {
    std::string text;

    // loop over all video files
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, 1920, 1080, GL_RGBA, GL_UNSIGNED_BYTE, image_.image.data());

    if (first_time_) {
        start_ = std::chrono::high_resolution_clock::now();
        last_ = start_;
        text = mmgraph_ptr->Convenience().SerializeAllParameters();
        old_param_text_ = text;
        first_time_ = false;
    } else {
        auto new_params = mmgraph_ptr->Convenience().SerializeAllParameters();
        text = parameter_diff(old_param_text_, new_params);
        old_param_text_ = new_params;
    }
    for (auto& [filename, stream_ctx] : stream_ctx_map_) {
        //auto& vid_ctx = stream_ctx_[0];
        auto& vid_ctx = stream_ctx[0];

        flipImage2RGB(vid_ctx, image_);

        rgb2yuv(vid_ctx);

        auto const current = std::chrono::high_resolution_clock::now();
        auto const time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(current - start_).count();
        //vid_ctx.yuvpic->pts = counter++;
        vid_ctx.yuvpic->pts = time_in_ms;

        if (!text.empty()) {
            auto base_ts = convert_to_timestamp(std::chrono::duration_cast<std::chrono::milliseconds>(last_ - start_));
            auto next_ts =
                convert_to_timestamp(std::chrono::duration_cast<std::chrono::milliseconds>(current - start_));

            //write_srt_entry(srt_file_, counter, base_ts, next_ts, std::string("Frame ") + std::to_string(counter++));
            write_srt_entry(srt_file_, counter++, base_ts, next_ts, text);
        }

        last_ = current;

        encodeFrame(vid_ctx);
    }
}


void megamol::frontend::Video_Service::fill_lua_callbacks() {
    frontend_resources::LuaCallbacksCollection callbacks;

    callbacks.add<frontend_resources::LuaCallbacksCollection::VoidResult, std::string>("mmStartVideoRecord",
        "(string filename)",
        {[&](std::string const& filename) -> frontend_resources::LuaCallbacksCollection::VoidResult {
            start_video_rec(filename);
            return frontend_resources::LuaCallbacksCollection::VoidResult{};
        }});

    callbacks.add<frontend_resources::LuaCallbacksCollection::VoidResult, std::string>("mmStopVideoRecord",
        "(string filename)",
        {[&](std::string const& filename) -> frontend_resources::LuaCallbacksCollection::VoidResult {
            // write subtitles into file
            // flush encoder and clean up
            stop_video_rec(filename);
            return frontend_resources::LuaCallbacksCollection::VoidResult{};
        }});

    // load video

    // pause video

    // unload video
}


void megamol::frontend::Video_Service::start_video_rec(std::string const& filename) {
    std::vector<StreamContext> sc;
    setup_video(filename, {1920, 1080}, sc);
    //setup_subtitles(filename, sc);
    stream_ctx_map_[filename] = std::move(sc);
}


void megamol::frontend::Video_Service::stop_video_rec(std::string const& filename) {
    auto fit = stream_ctx_map_.find(filename);
    if (fit != stream_ctx_map_.end()) {
        // write subtitles
        //encode_sub(fit->second);
        flush_encoder(fit->second[0]);

        srt_file_.close();
        setup_subtitles("./test_out.srt", fit->second);

        encode_sub(fit->second);

        auto ret = av_write_trailer(fit->second[0].fmt_ctx);
        // flush encoder
        // clean up
    }
}
