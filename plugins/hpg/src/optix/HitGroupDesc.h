#pragma once

#include <memory>

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"

#include "hpg/optix/CallProgramGroupDesc.h"

#include "optix.h"

namespace megamol::hpg::optix {

class HitGroupDesc : public core::Module {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() { return "HitGroupDesc"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() { return "Module representing an OptiX hit group descriptor"; }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() { return true; }

    HitGroupDesc();

    virtual ~HitGroupDesc();

protected:
    bool create() override;

    void release() override;

private:
    bool get_desc_cb(core::Call& c);

    core::CalleeSlot _out_desc_slot;

    core::CallerSlot _in_shading_slot;

    core::CallerSlot _in_intersection_slot;

    std::shared_ptr<OptixProgramGroupDesc> _desc;
};

} // end namespace megamol::hpg::optix