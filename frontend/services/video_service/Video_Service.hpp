#pragma once

#include <chrono>
#include <fstream>
#include <unordered_map>

#include "AbstractFrontendService.hpp"

#include "Screenshots.h"
#include "mmcore/MegaMolGraph.h"
#include "GUIRegisterWindow.h"
#include "ImageWrapper.h"

#include "glowl/Texture2D.hpp"

namespace megamol::frontend {

struct StreamContext;

class Video_Service : public AbstractFrontendService {
public:
    std::string serviceName() const override {
        return "Video_Service";
    }

    Video_Service();
    ~Video_Service();

    bool init(void* configPtr) override;
    void close() override;

    std::vector<FrontendResource>& getProvidedResources() override;

    const std::vector<std::string> getRequestedResourceNames() const override;
    void setRequestedResources(std::vector<FrontendResource> resources) override;

    void updateProvidedResources() override;

    void digestChangedRequestedResources() override;

    void resetProvidedResources() override;

    // check whether recording is started or stopped
    void preGraphRender() override;

    // encode current frame
    void postGraphRender() override;

protected:
private:
    using LuaFuncType = std::function<std::tuple<bool, std::string>(std::string const&)>;

    void fill_lua_callbacks();

    void start_video_rec(std::string const& filename);

    void stop_video_rec(std::string const& filename);

    void create_playback_window(megamol::frontend_resources::ImageWrapper const& image);

    std::vector<std::string> requestedResourcesNames_;

    std::vector<megamol::frontend::FrontendResource> providedResources_;

    megamol::frontend_resources::ScreenshotImageData image_;

    int counter = 0;

    std::unordered_map<std::string, std::vector<StreamContext>> stream_ctx_map_;

    bool first_time_ = true;

    std::chrono::high_resolution_clock::time_point start_;

    std::chrono::high_resolution_clock::time_point last_;

    std::ofstream srt_file_;

    megamol::core::MegaMolGraph* mmgraph_ptr = nullptr;

    megamol::frontend_resources::GUIRegisterWindow* guireg_ptr = nullptr;

    //std::shared_ptr<glowl::Texture2D> ogl_texture_;

    std::shared_ptr<frontend_resources::ImageWrapper> iw_;

    GLuint ogl_texture_ = 0;

    LuaFuncType const* execute_lua_;

    std::string old_param_text_;
};
} // namespace megamol::frontend
