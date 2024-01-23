#pragma once

#include <vector>

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"

#include "mmcore/param/ParamSlot.h"

#include "geometry_calls/MultiParticleDataCall.h"

#include "optix/Context.h"
#include "optix/Utils.h"

#include <cuda.h>

namespace megamol::optix_hpg {
class PKDGeometry : public core::Module {
public:
    enum class PKDMode { STANDARD, TREELETS };

    static const char* ClassName(void) {
        return "PKDGeometry";
    }

    static const char* Description(void) {
        return "PKD Geometry for OptiX";
    }

    static bool IsAvailable(void) {
        return true;
    }

    PKDGeometry();

    virtual ~PKDGeometry();

protected:
    bool create() override;

    void release() override;

private:
    bool get_data_cb(core::Call& c);

    bool get_extents_cb(core::Call& c);

    bool init(Context const& ctx);

    bool assert_data(geocalls::MultiParticleDataCall const& call, Context const& ctx);

    void free_data(Context const& ctx) {
        for (auto& el : particle_data_) {
            CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
        }
        for (auto& el : radius_data_) {
            CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
        }
        for (auto& el : color_data_) {
            CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
        }
    }

    core::CalleeSlot out_geo_slot_;

    core::CallerSlot in_data_slot_;

    core::param::ParamSlot mode_slot_;

    core::param::ParamSlot compression_slot_;

    core::param::ParamSlot threshold_slot_;

    std::vector<CUdeviceptr> particle_data_;

    std::vector<CUdeviceptr> radius_data_;

    std::vector<CUdeviceptr> color_data_;
};
} // namespace megamol::optix_hpg
