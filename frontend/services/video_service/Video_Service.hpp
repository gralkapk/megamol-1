#pragma once

#include "AbstractFrontendService.hpp"

namespace megamol::frontend {
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
};
} // namespace megamol::frontend
