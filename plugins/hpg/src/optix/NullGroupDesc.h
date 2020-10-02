#pragma once

#include <memory>

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"

#include "hpg/optix/CallProgramGroupDesc.h"

#include "optix.h"

namespace megamol::hpg::optix {

class NullGroupDesc : public core::Module {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() { return "NullGroupDesc"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() { return "Module providing empty program group descriptor"; }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() { return true; }

    NullGroupDesc();

    virtual ~NullGroupDesc();

protected:
    bool create() override;

    void release() override;

private:
    bool get_desc_cb(core::Call& c);

    core::CalleeSlot _out_desc_slot;

    std::shared_ptr<OptixProgramGroupDesc> _desc;
}; // end class NullGroupDesc

} // end namespace megamol::hpg::optix
