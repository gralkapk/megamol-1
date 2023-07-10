#include "Video_Service.hpp"

#include <regex>
#include <sstream>
#include <filesystem>

#include "video_util.hpp"

#include "glad/gl.h"

#include <imgui_stdlib.h>


#include "ImageWrapper.h"
#include "LuaCallbacksCollection.h"

constexpr bool writeVideo = true;

megamol::frontend::Video_Service::Video_Service() {}


megamol::frontend::Video_Service::~Video_Service() {}


bool megamol::frontend::Video_Service::init(void* configPtr) {

    //requestedResourcesNames_ = {"RegisterLuaCallback", "MegaMolGraph"};
    requestedResourcesNames_ = {
        "MegaMolGraph", "GUIRegisterWindow", "ExecuteLuaScript", "SetScriptPath", "FramebufferEvents"};


    // for test purposes
    

    //image_.resize(1920, 1080);

    //glGenTextures(1, &ogl_texture_);
    //glBindTexture(GL_TEXTURE_2D, ogl_texture_);

    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, 1920, 1080, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);

    ///*ogl_texture_ = std::make_shared<glowl::Texture2D>(
    //    "playback_texture", glowl::TextureLayout(GL_RGB8, 1920, 1080, 0, GL_RGB, GL_UNSIGNED_BYTE, 1), nullptr);*/

    //iw_ = std::make_shared<frontend_resources::ImageWrapper>(frontend_resources::wrap_image(
    //    {1920, 1080}, ogl_texture_, frontend_resources::ImageWrapper::DataChannels::RGB8));

    resize();

    return true;
}


void megamol::frontend::Video_Service::resize() {
    image_.resize(fbo_size_.x, fbo_size_.y);

    if (glIsTexture(ogl_texture_)) {
        glDeleteTextures(1, &ogl_texture_);
    }

    glGenTextures(1, &ogl_texture_);
    glBindTexture(GL_TEXTURE_2D, ogl_texture_);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, fbo_size_.x, fbo_size_.y, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);

    /*ogl_texture_ = std::make_shared<glowl::Texture2D>(
        "playback_texture", glowl::TextureLayout(GL_RGB8, 1920, 1080, 0, GL_RGB, GL_UNSIGNED_BYTE, 1), nullptr);*/

    iw_ = std::make_shared<frontend_resources::ImageWrapper>(frontend_resources::wrap_image(
        frontend_resources::ImageWrapper::ImageSize{static_cast<size_t>(fbo_size_.x), static_cast<size_t>(fbo_size_.y)},
        ogl_texture_, frontend_resources::ImageWrapper::DataChannels::RGB8));
}


void megamol::frontend::Video_Service::close() {
    if constexpr (writeVideo) {
        stop_video_rec("./test_out.mkv");
    }

    glDeleteTextures(1, &ogl_texture_);
}


std::vector<megamol::frontend::FrontendResource>& megamol::frontend::Video_Service::getProvidedResources() {
    return providedResources_;
}


const std::vector<std::string> megamol::frontend::Video_Service::getRequestedResourceNames() const {
    return requestedResourcesNames_;
}


void megamol::frontend::Video_Service::setRequestedResources(std::vector<FrontendResource> resources) {
    mmgraph_ptr = const_cast<megamol::core::MegaMolGraph*>(&resources[0].getResource<megamol::core::MegaMolGraph>());
    guireg_ptr = const_cast<megamol::frontend_resources::GUIRegisterWindow*>(
        &resources[1].getResource<megamol::frontend_resources::GUIRegisterWindow>());
    execute_lua_ = &resources[2].getResource<LuaFuncType>();
    set_script_path_ = &resources[3].getResource<SetScriptPath>();
    fbo_events_ = &resources[4].getResource<megamol::frontend_resources::FramebufferEvents>();

    fbo_size_ = glm::ivec2(fbo_events_->previous_state.width, fbo_events_->previous_state.height);

    

    //if (fbo_size_ != glm::ivec2(1)) {
    //    if constexpr (writeVideo) {
    //        start_video_rec("./test_out.mkv");
    //        srt_file_ = std::ofstream("./test_out.srt");
    //    } else {
    //        std::vector<StreamContext> sc;
    //        open_video("./test_out.mkv", sc);
    //        stream_ctx_map_["./test_out.mkv"] = std::move(sc);
    //    }
    //}

    create_playback_window(*iw_.get());
    create_recorder_window();
}


void megamol::frontend::Video_Service::updateProvidedResources() {}


void megamol::frontend::Video_Service::digestChangedRequestedResources() {
    // check for resize
    if (fbo_events_->is_resized()) {
        auto current_size = fbo_events_->size_events.back();
        fbo_size_ = glm::ivec2(current_size.width, current_size.height);

        resize();

        if constexpr (writeVideo) {
            start_video_rec("./test_out.mkv");
            srt_file_ = std::ofstream("./test_out.srt");
        } else {
            std::vector<StreamContext> sc;
            open_video("./test_out.mkv", sc);
            stream_ctx_map_["./test_out.mkv"] = std::move(sc);
        }
    }
}


void megamol::frontend::Video_Service::resetProvidedResources() {}


void megamol::frontend::Video_Service::preGraphRender() {
    // handle start and stop
    // need lua callbacks for that
    // read frame
    if constexpr (!writeVideo) {
        auto& stream_ctx = stream_ctx_map_["./test_out.mkv"];
        auto& vid_ctx = stream_ctx[0];
        std::string test_txt;
        decode_frame(vid_ctx.in_fmt_ctx, stream_ctx, image_, test_txt);

        // write texture
        glBindTexture(GL_TEXTURE_2D, ogl_texture_);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, fbo_size_.x, fbo_size_.y, GL_RGBA, GL_UNSIGNED_BYTE, image_.image.data());
        glBindTexture(GL_TEXTURE_2D, 0);
        if (!test_txt.empty()) {
            // 0,0,Default,,0,0,0,,
            test_txt = std::regex_replace(test_txt, std::regex("\\\\N"), "\n");
            test_txt = std::regex_replace(test_txt, std::regex("\\d+.\\d.Default..\\d.\\d.\\d.."), "");
            std::cout << "[SRT]\n" << test_txt << "\n[SRT]" << std::endl;
            //set_script_path_->operator()("C:\\data\\dev\\vidmol\\build\\x64-Debug\\install\\examples\\testspheres.lua");
            auto result = execute_lua_->operator()(test_txt);
            if (!std::get<0>(result)) {
                std::cout << "[LUA ERROR] " << std::get<1>(result) << std::endl;
            }
            if (test_txt.find("mmCreateView") != std::string::npos) {
                std::cout << "[SRT] Found View" << std::endl;
                //execute_lua_->operator()("mmRenderNextFrame()");
            }
        }
    }
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
    if constexpr (writeVideo) {
        std::string text;

        // loop over all video files
        glReadBuffer(GL_FRONT);
        glReadPixels(0, 0, fbo_size_.x, fbo_size_.y, GL_RGBA, GL_UNSIGNED_BYTE, image_.image.data());

        for (auto& [filename, vc] : video_ctx_map_) {
            capture_frame(vc, image_);
        }

        if (first_time_) {
            start_ = std::chrono::high_resolution_clock::now();
            last_ = start_;
            //text = mmgraph_ptr->Convenience().SerializeAllParameters();
            text = mmgraph_ptr->Convenience().SerializeGraph();
            old_param_text_ = text;
            first_time_ = false;
        } else {
            //auto new_params = mmgraph_ptr->Convenience().SerializeAllParameters();
            auto new_params = mmgraph_ptr->Convenience().SerializeGraph();
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
                auto base_ts =
                    convert_to_timestamp(std::chrono::duration_cast<std::chrono::milliseconds>(last_ - start_));
                auto next_ts =
                    convert_to_timestamp(std::chrono::duration_cast<std::chrono::milliseconds>(current - start_));

                //write_srt_entry(srt_file_, counter, base_ts, next_ts, std::string("Frame ") + std::to_string(counter++));
                write_srt_entry(srt_file_, counter++, base_ts, next_ts, text);
            }

            last_ = current;

            encodeFrame(vid_ctx);
        }
    } else {
#if 0
        // read frame
        auto& stream_ctx = stream_ctx_map_["./test_out.mkv"];
        auto& vid_ctx = stream_ctx[0];
        std::string test_txt;
        decode_frame(vid_ctx.in_fmt_ctx, stream_ctx, image_, test_txt);

        // write texture
        glBindTexture(GL_TEXTURE_2D, ogl_texture_);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 1920, 1080, GL_RGBA, GL_UNSIGNED_BYTE, image_.image.data());
        glBindTexture(GL_TEXTURE_2D, 0);
        if (!test_txt.empty()) {
            // 0,0,Default,,0,0,0,,
            test_txt = std::regex_replace(test_txt, std::regex("\\\\N"), "\n");
            test_txt = std::regex_replace(test_txt, std::regex("\\d.\\d.Default..\\d.\\d.\\d.."), "");
            std::cout << "[SRT]\n" << test_txt << "\n[SRT]" << std::endl;
            set_script_path_->operator()("C:\\data\\dev\\vidmol\\build\\x64-Debug\\install\\examples\\testspheres.lua");
            auto result = execute_lua_->operator()(test_txt);
            if (!std::get<0>(result)) {
                std::cout << "[LUA ERROR] " << std::get<1>(result) << std::endl; 
            }
            if (test_txt.find("mmCreateView") != std::string::npos) {
                std::cout << "[SRT] Found View" << std::endl;
                execute_lua_->operator()("mmRenderNextFrame()");
            }
        }
#endif
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

    // register sink with video output
}


void megamol::frontend::Video_Service::start_video_rec(std::string const& filename) {
    std::vector<StreamContext> sc;
    setup_video(filename, fbo_size_, sc);

    VideoContext vc;
    vc.dim = fbo_size_;
    vc.stream_ctx = sc;
    std::filesystem::path p(filename);
    auto const full_p = std::filesystem::temp_directory_path() / p.filename().replace_extension("srt");
    vc.srt_file_path = full_p;
    vc.srt_file.open(vc.srt_file_path);
    video_ctx_map_[filename] = std::move(vc);

    //stream_ctx_map_[filename] = std::move(sc);
}


void megamol::frontend::Video_Service::stop_video_rec(std::string const& filename) {
    //auto fit = stream_ctx_map_.find(filename);
    //if (fit != stream_ctx_map_.end()) {
    //    // write subtitles
    //    //encode_sub(fit->second);
    //    flush_encoder(fit->second[0]);

    //    srt_file_.close();
    //    setup_subtitles("./test_out.srt", fit->second);

    //    encode_sub(fit->second);

    //    auto ret = av_write_trailer(fit->second[0].fmt_ctx);
    //    // flush encoder
    //    // clean up
    //}

    auto v_fit = video_ctx_map_.find(filename);
    if (v_fit != video_ctx_map_.end()) {
        flush_encoder(v_fit->second.stream_ctx[0]);

        v_fit->second.srt_file.close();
        setup_subtitles(v_fit->second.srt_file_path.string(), v_fit->second.stream_ctx);

        encode_sub(v_fit->second.stream_ctx);

        auto ret = av_write_trailer(v_fit->second.stream_ctx[0].fmt_ctx);
    }

    video_ctx_map_.erase(filename);
}


struct iw_functor {
    void init(std::vector<megamol::frontend_resources::ImageWrapper> const& images) {
        images_ = images;
    }
    std::vector<megamol::frontend_resources::ImageWrapper> images_;
};


void megamol::frontend::Video_Service::create_playback_window(megamol::frontend_resources::ImageWrapper const& image) {
    auto iw = std::make_shared<iw_functor>();
    iw->init({image});

    auto win_func = std::bind(
        [](std::shared_ptr<iw_functor> iw, megamol::gui::AbstractWindow::BasicConfig& window_config) {
            window_config.flags = ImGuiWindowFlags_AlwaysAutoResize;

            for (auto& image : iw->images_) {
                ImGui::Image(image.referenced_image_handle, ImVec2{(float)image.size.width, (float)image.size.height},
                    ImVec2(0, 1), ImVec2(1, 0));
            }
        },
        iw, std::placeholders::_1);

    guireg_ptr->register_window("Video: playback", win_func);
}


void megamol::frontend::Video_Service::create_recorder_window() {
    auto win_func = [&](megamol::gui::AbstractWindow::BasicConfig& window_config) {
        // record button
        // stop button
        // output file
        // capture entrypoint
        std::string output_file;
        ImGui::InputText("ouput", &output_file);

        if (ImGui::Button("record")) {
            // start record
            if (!output_file.empty()) {
                start_video_rec(output_file);
            } else {
                core::utility::log::Log::DefaultLog.WriteError("[Video_Service] Provide path to output file first");
            }
        }
        if (ImGui::Button("stop")) {
            // stop record
            if (!output_file.empty()) {
                stop_video_rec(output_file);
            } else {
                core::utility::log::Log::DefaultLog.WriteError("[Video_Service] Provide path to output file first");
            }
        }
    };

    guireg_ptr->register_window("Video: recorder", win_func);
}


void megamol::frontend::Video_Service::capture_frame(
    VideoContext& vc, megamol::frontend_resources::ScreenshotImageData const& image) {
    std::string text;
    if (vc.first_time) {
        vc.start = std::chrono::high_resolution_clock::now();
        vc.last = vc.start;
        //text = mmgraph_ptr->Convenience().SerializeAllParameters();
        text = mmgraph_ptr->Convenience().SerializeGraph();
        vc.old_param_list = text;
        vc.first_time = false;
    } else {
        //auto new_params = mmgraph_ptr->Convenience().SerializeAllParameters();
        auto new_params = mmgraph_ptr->Convenience().SerializeGraph();
        text = parameter_diff(vc.old_param_list, new_params);
        vc.old_param_list = new_params;
    }

    //auto& vid_ctx = stream_ctx_[0];
    auto& vid_ctx = vc.stream_ctx[0];

    flipImage2RGB(vid_ctx, image);

    rgb2yuv(vid_ctx);

    auto const current = std::chrono::high_resolution_clock::now();
    auto const time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(current - vc.start).count();
    //vid_ctx.yuvpic->pts = counter++;
    vid_ctx.yuvpic->pts = time_in_ms;

    if (!text.empty()) {
        auto base_ts = convert_to_timestamp(std::chrono::duration_cast<std::chrono::milliseconds>(vc.last - vc.start));
        auto next_ts = convert_to_timestamp(std::chrono::duration_cast<std::chrono::milliseconds>(current - vc.start));

        //write_srt_entry(srt_file_, counter, base_ts, next_ts, std::string("Frame ") + std::to_string(counter++));
        write_srt_entry(vc.srt_file, vc.counter++, base_ts, next_ts, text);
    }

    vc.last = current;

    encodeFrame(vid_ctx);
}
