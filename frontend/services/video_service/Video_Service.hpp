#pragma once

#include <unordered_map>

#include "AbstractFrontendService.hpp"

#include "Screenshots.h"

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
    void fill_lua_callbacks();

    void start_video_rec(std::string const& filename);

    void stop_video_rec(std::string const& filename);

    std::vector<std::string> requestedResourcesNames_;

    std::vector<megamol::frontend::FrontendResource> providedResources_;

    megamol::frontend_resources::ScreenshotImageData image_;

    int counter = 0;

    std::unordered_map<std::string, std::vector<StreamContext>> stream_ctx_map_;
};
} // namespace megamol::frontend
