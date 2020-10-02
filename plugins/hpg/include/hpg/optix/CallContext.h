#pragma once

#include <memory>

#include "mmcore/Call.h"
#include "mmcore/factories/CallAutoDescription.h"

#include "optix.h"

namespace megamol::hpg::optix {

class CallContext : public core::Call {
public:
    static const char* ClassName(void) { return "CallContext"; }

    static const char* Description(void) { return "Transports an OptiX context"; }

    static unsigned int FunctionCount() { return 1; }

    static const char* FunctionName(unsigned int idx) {
        switch (idx) {
        case 0:
            return "GetContext";
        }
        return nullptr;
    }

    CallContext() = default;

    virtual ~CallContext() = default;

    bool operator()(unsigned int func = 0) = delete;

    bool is_dirty() { return _dirty; }

    void set_dirty() { _dirty = true; }

    void reset_dirty() { _dirty = false; }

    void set_context(std::shared_ptr<OptixDeviceContext> const& ctx) { _ctx = ctx; }

    std::shared_ptr<OptixDeviceContext> get_context() {
        if (!dynamic_cast<Call*>(this)->operator()(0)) {
            return nullptr;
        }

        return _ctx;
    }

private:
    bool _dirty;

    std::shared_ptr<OptixDeviceContext> _ctx;

}; // end class CallContext

using CallContextDescription = megamol::core::factories::CallAutoDescription<CallContext>;

} // end namespace megamol::hpg::optix
