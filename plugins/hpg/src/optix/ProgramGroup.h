#pragma once

#include <memory>
#include <vector>

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"

#include "hpg/optix/CallProgramGroupDesc.h"

#include "optix.h"

namespace megamol::hpg::optix {

class ProgramGroup : public core::Module {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() { return "ProgramGroup"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() { return "Module representing an OptiX program group"; }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() { return true; }

    ProgramGroup();

    virtual ~ProgramGroup();

protected:
    bool create() override;

    void release() override;

private:
    bool get_desc_cb(core::Call& c);

    core::CalleeSlot _out_groups_slot;

    core::CallerSlot _in_context_slot;

    core::CallerSlot _in_raygen_slot;

    core::CallerSlot _in_miss_slot;

    core::CallerSlot _in_exception_slot;

    core::CallerSlot _in_hitgroup_slot;

    core::CallerSlot _in_callables_slot;

    std::shared_ptr<std::vector<OptixProgramGroup>> _groups;
}; // end class ProgramGroup

} // end namespace megamol::hpg::optix