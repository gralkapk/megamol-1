#pragma once

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"

#include "geometry_calls/MultiParticleDataCall.h"

namespace megamol::moldyn {
class Accumulator : public core::Module {
public:
    Accumulator();

    virtual ~Accumulator();

protected:
    bool create() override;

    void release() override;

private:
    enum class dir_t {
        BACKWARD,
        CENTRAL,
        FORWARD
    };

    bool get_data_cb(core::Call& c);

    bool get_extent_cb(core::Call& c);

    bool isDirty() {
        return window_size_slot_.IsDirty();
    }

    void resetDirty() {
        window_size_slot_.ResetDirty();
    }

    void collect_backward(geocalls::MultiParticleDataCall* data, int window_size);

    void collect_central(geocalls::MultiParticleDataCall* data, int window_size);

    void collect_forward(geocalls::MultiParticleDataCall* data, int window_size);

    core::CalleeSlot data_out_slot_;

    core::CallerSlot data_in_slot_;

    core::param::ParamSlot window_size_slot_;

    core::param::ParamSlot direction_slot_;

    unsigned int frame_id_ = 0;

    uint64_t in_data_hash_ = std::numeric_limits<uint64_t>::max();

    uint64_t out_data_hash_ = 0;
};
} // namespace megamol::moldyn
