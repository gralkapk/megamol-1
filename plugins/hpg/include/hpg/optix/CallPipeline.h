#pragma once

#include <memory>

#include "mmcore/Call.h"
#include "mmcore/factories/CallAutoDescription.h"

#include "optix.h"

namespace megamol::hpg::optix {

class CallPipeline : public core::Call {
public:
    static const char* ClassName(void) { return "CallPipeline"; }

    static const char* Description(void) { return "Transports an OptiX pipeline"; }

    static unsigned int FunctionCount() { return 1; }

    static const char* FunctionName(unsigned int idx) {
        switch (idx) {
        case 0:
            return "GetPipeline";
        }
        return nullptr;
    }

    CallPipeline() = default;

    virtual ~CallPipeline() = default;

    bool operator()(unsigned int func = 0) = delete;

    bool is_dirty() { return _dirty; }

    void set_dirty() { _dirty = true; }

    void reset_dirty() { _dirty = false; }

    void set_pipeline(std::shared_ptr<OptixPipeline> const& pipe) { _pipe = pipe; }

    std::shared_ptr<OptixPipeline> get_pipeline() {
        if (!dynamic_cast<Call*>(this)->operator()(0)) {
            return nullptr;
        }

        return _pipe;
    }

private:
    bool _dirty;

    std::shared_ptr<OptixPipeline> _pipe;

}; // end class CallPipeline

using CallPipelineDescription = megamol::core::factories::CallAutoDescription<CallPipeline>;

} // end namespace megamol::hpg::optix
