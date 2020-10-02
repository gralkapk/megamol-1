#pragma once

#include <memory>

#include "mmcore/Call.h"

#include "optix.h"

namespace megamol::hpg::optix {

class CallProgramGroupDesc : public core::Call {
public:
    static const char* ClassName(void) { return "CallProgramGroupDesc"; }

    static const char* Description(void) { return "Transports an OptiX program group descriptor"; }

    static unsigned int FunctionCount() { return 1; }

    static const char* FunctionName(unsigned int idx) {
        switch (idx) {
        case 0:
            return "GetDescriptor";
        }
        return nullptr;
    }

    CallProgramGroupDesc() = default;

    virtual ~CallProgramGroupDesc() = default;

    bool operator()(unsigned int func = 0) = delete;

    bool is_dirty() {
        return _dirty;
    }

    void set_dirty() { _dirty = true; }

    void reset_dirty() { _dirty = false; }

    void set_descriptor(std::shared_ptr<OptixProgramGroupDesc> const& desc) { _desc = desc; }

    std::shared_ptr<OptixProgramGroupDesc const> get_descriptor() {
        if (!dynamic_cast<Call*>(this)->operator()(0)) {
            return nullptr;
        }

        return _desc;
    }

private:
    bool _dirty;

    std::shared_ptr<OptixProgramGroupDesc const> _desc;

}; // end class CallProgramGroupDesc

using CallProgramGroupDescDescription = megamol::core::factories::CallAutoDescription<CallProgramGroupDesc>;

} // end namespace megamol::hpg::optix
