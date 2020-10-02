#pragma once

#include <memory>
#include <vector>

#include "mmcore/Call.h"
#include "mmcore/factories/CallAutoDescription.h"

#include "optix.h"

namespace megamol::hpg::optix {

class CallProgramGroup : public core::Call {
public:
    static const char* ClassName(void) { return "CallProgramGroup"; }

    static const char* Description(void) { return "Transports OptiX program groups"; }

    static unsigned int FunctionCount() { return 1; }

    static const char* FunctionName(unsigned int idx) {
        switch (idx) {
        case 0:
            return "GetContext";
        }
        return nullptr;
    }

    CallProgramGroup() = default;

    virtual ~CallProgramGroup() = default;

    bool operator()(unsigned int func = 0) = delete;

    bool is_dirty() { return _dirty; }

    void set_dirty() { _dirty = true; }

    void reset_dirty() { _dirty = false; }

    void set_groups(std::shared_ptr<std::vector<OptixProgramGroup>> const& groups) { _groups = groups; }

    std::shared_ptr<std::vector<OptixProgramGroup>> get_groups() {
        if (!dynamic_cast<Call*>(this)->operator()(0)) {
            return nullptr;
        }

        return _groups;
    }

private:
    bool _dirty;

    std::shared_ptr<std::vector<OptixProgramGroup>> _groups;

}; // end class CallProgramGroup

using CallProgramGroupDescription = megamol::core::factories::CallAutoDescription<CallProgramGroup>;

} // end namespace megamol::hpg::optix
