#pragma once

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"

namespace megamol::thermodyn {
class PCLSurface : public core::Module {
public:
    /** Return module class name */
    static const char* ClassName(void) {
        return "PCLSurface";
    }

    /** Return module class description */
    static const char* Description(void) {
        return "PCLSurface";
    }

    /** Module is always available */
    static bool IsAvailable(void) {
        return true;
    }

    PCLSurface();

    virtual ~PCLSurface();

protected:
    bool create() override;

    void release() override;

private:
    bool get_data_cb(core::Call& c);

    bool get_extent_cb(core::Call& c);

    core::CalleeSlot out_data_slot_;

    core::CallerSlot in_data_slot_;

    uint64_t in_data_hash_ = std::numeric_limits<uint64_t>::max();

    uint64_t out_data_hash_ = 0;

    unsigned int frame_id_ = 0;
};
} // namespace megamol::thermodyn
