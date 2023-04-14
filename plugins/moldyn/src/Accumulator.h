#pragma once

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"

#include "geometry_calls/MultiParticleDataCall.h"

#include <tuple>

#include "glm/glm.hpp"

namespace megamol::moldyn {
class Accumulator : public core::Module {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() {
        return "Accumulator";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() {
        return "Accumulator";
    }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() {
        return true;
    }

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
        return window_size_slot_.IsDirty() || direction_slot_.IsDirty();
    }

    void resetDirty() {
        window_size_slot_.ResetDirty();
        direction_slot_.ResetDirty();
    }

    std::tuple<std::vector<std::vector<uint64_t>>, std::vector<std::vector<glm::vec3>>,
        std::vector<std::vector<glm::vec3>>, std::vector<std::vector<glm::vec4>>>
    collect_backward(geocalls::MultiParticleDataCall* data, int window_size);

    void collect_central(geocalls::MultiParticleDataCall* data, int window_size);

    void collect_forward(geocalls::MultiParticleDataCall* data, int window_size);

    core::CalleeSlot data_out_slot_;

    core::CallerSlot data_in_slot_;

    core::param::ParamSlot window_size_slot_;

    core::param::ParamSlot direction_slot_;

    unsigned int frame_id_ = 0;

    uint64_t in_data_hash_ = std::numeric_limits<uint64_t>::max();

    uint64_t out_data_hash_ = 0;

    std::vector<std::vector<glm::vec3>> avg_pos_;
    std::vector<std::vector<glm::vec3>> avg_dir_;
    std::vector<std::vector<glm::vec4>> avg_col_;
    std::vector<std::vector<uint64_t>> base_id_;
    std::vector<float> global_radii_;
};
} // namespace megamol::moldyn
