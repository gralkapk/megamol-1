#include "Video_Service.hpp"

#include "video_util.hpp"

#include "glad/gl.h"


#include "LuaCallbacksCollection.h"


megamol::frontend::Video_Service::Video_Service() {}


megamol::frontend::Video_Service::~Video_Service() {}


bool megamol::frontend::Video_Service::init(void* configPtr) {

    requestedResourcesNames_ = {"RegisterLuaCallback"};

    return true;
}


void megamol::frontend::Video_Service::close() {}


std::vector<megamol::frontend::FrontendResource>& megamol::frontend::Video_Service::getProvidedResources() {
    // TODO: insert return statement here
}


const std::vector<std::string> megamol::frontend::Video_Service::getRequestedResourceNames() const {
    return std::vector<std::string>();
}


void megamol::frontend::Video_Service::setRequestedResources(std::vector<FrontendResource> resources) {}


void megamol::frontend::Video_Service::updateProvidedResources() {}


void megamol::frontend::Video_Service::digestChangedRequestedResources() {}


void megamol::frontend::Video_Service::resetProvidedResources() {}


void megamol::frontend::Video_Service::preGraphRender() {
    // handle start and stop
    // need lua callbacks for that
}


void megamol::frontend::Video_Service::postGraphRender() {
    // loop over all video files

    auto& vid_ctx = stream_ctx_[0];

    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, vid_ctx.dim.x, vid_ctx.dim.y, GL_RGBA, GL_UNSIGNED_BYTE, image_.image.data());

    flipRGB(vid_ctx, image_);

    rgb2yuv(vid_ctx);

    vid_ctx.yuvpic->pts = counter++;

    //encodeFrame(vid_ctx)
}


void megamol::frontend::Video_Service::fill_lua_callbacks() {
    frontend_resources::LuaCallbacksCollection callbacks;

    callbacks.add<frontend_resources::LuaCallbacksCollection::VoidResult, std::string>("mmStartVideoRecord",
        "(string filename)",
        {[&](std::string const& filename) -> frontend_resources::LuaCallbacksCollection::VoidResult {
            std::vector<StreamContext> sc;
            setup_video(filename, {1920, 1080}, sc);
            stream_ctx_map_[filename] = sc;
            return frontend_resources::LuaCallbacksCollection::VoidResult{};
        }});

    callbacks.add<frontend_resources::LuaCallbacksCollection::VoidResult, std::string>("mmStopVideoRecord",
        "(string filename)",
        {[](std::string const& filename) -> frontend_resources::LuaCallbacksCollection::VoidResult {
            // write subtitles into file
            // flush encoder and clean up
            return frontend_resources::LuaCallbacksCollection::VoidResult{};
        }});
}
