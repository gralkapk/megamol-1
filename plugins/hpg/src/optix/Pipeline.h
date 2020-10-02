#pragma once

#include <memory>

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"

#include "optix.h"

namespace megamol::hpg::optix {

class Pipeline : public core::Module {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() { return "Pipeline"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() { return "Module representing an OptiX pipeline"; }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() { return true; }

    Pipeline();

    virtual ~Pipeline();

protected:
    bool create() override;

    void release() override;

private:
    bool get_pipeline_cb(core::Call& c);

    core::CalleeSlot _out_pipeline_slot;

    core::CallerSlot _in_context_slot;

    core::CallerSlot _in_groups_slot;

    std::shared_ptr<OptixPipeline> _pipe;
}; // end class Pipeline

} // end namespace megamol::hpg::optix