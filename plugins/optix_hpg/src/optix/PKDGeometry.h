// PKD implementation:
// ======================================================================== //
// Copyright 2018-2019 Ingo Wald                                            //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

// ======================================================================== //
// Modified 2019-2022 VISUS - University of Stuttgart                       //
// ======================================================================== //

#pragma once

#include <vector>

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"

#include "mmcore/param/ParamSlot.h"

#include "geometry_calls/MultiParticleDataCall.h"

#include "SBT.h"
#include "optix/Context.h"
#include "optix/Utils.h"
#include "MMOptixModule.h"

#include <cuda.h>

#include "pkd.h"

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

    bool createSBTRecords(geocalls::MultiParticleDataCall const& call, Context const& ctx);

    bool has_color(geocalls::SimpleSphericalParticles const& parts) const {
        auto const color_type = parts.GetColourDataType();
        return (color_type != geocalls::SimpleSphericalParticles::COLDATA_NONE) &&
               (color_type != geocalls::SimpleSphericalParticles::COLDATA_DOUBLE_I) &&
               (color_type != geocalls::SimpleSphericalParticles::COLDATA_FLOAT_I);
    }

    bool has_global_radius(geocalls::SimpleSphericalParticles const& parts) const {
        auto const vert_type = parts.GetVertexDataType();
        return vert_type != geocalls::SimpleSphericalParticles::VERTDATA_FLOAT_XYZR;
    }

    void free_data(Context const& ctx) {
        for (auto& el : particle_data_) {
            CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
        }
        /*for (auto& el : radius_data_) {
            CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
        }*/
        for (auto& el : color_data_) {
            CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
        }
        for (auto& el : treelets_data_) {
            CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
        }
    }

    core::CalleeSlot out_geo_slot_;

    core::CallerSlot in_data_slot_;

    core::param::ParamSlot mode_slot_;

    core::param::ParamSlot compression_slot_;

    core::param::ParamSlot threshold_slot_;

    std::vector<SBTRecord<device::PKDGeoData>> sbt_records_;

    std::vector<SBTRecord<device::TreeletsGeoData>> treelets_sbt_records_;

    std::array<OptixProgramGroup, 2> program_groups_;

    std::vector<CUdeviceptr> particle_data_;

    //std::vector<CUdeviceptr> radius_data_;

    std::vector<CUdeviceptr> color_data_;

    std::vector<CUdeviceptr> treelets_data_;

    std::vector<box3f> local_boxes_;

    CUdeviceptr geo_buffer_ = 0;

    OptixTraversableHandle geo_handle_;

    MMOptixModule pkd_module_;

    MMOptixModule pkd_occlusion_module_;

    MMOptixModule treelets_module_;

    MMOptixModule treelets_occlusion_module_;

    uint64_t sbt_version = 0;

    uint64_t program_version = 0;

    uint64_t geo_version = 0;

    unsigned int frame_id_ = std::numeric_limits<unsigned int>::max();

    std::size_t data_hash_ = std::numeric_limits<std::size_t>::max();
};

inline size_t parent(size_t C) {
    if (C == 0)
        return 0;
    return (C - 1) / 2;
}

inline size_t lChild(size_t P) {
    return 2 * P + 1;
}
inline size_t rChild(size_t P) {
    return 2 * P + 2;
}

template<class Comp>
inline void trickle(const Comp& worse, size_t P, device::PKDParticle* particle, size_t N, int dim) {
    if (P >= N)
        return;

    while (1) {
        const size_t L = lChild(P);
        const size_t R = rChild(P);
        const bool lValid = (L < N);
        const bool rValid = (R < N);

        if (!lValid)
            return;
        size_t C = L;
        if (rValid && worse(particle[R].pos[dim], particle[L].pos[dim]))
            C = R;

        if (!worse(particle[C].pos[dim], particle[P].pos[dim]))
            return;

        std::swap(particle[C], particle[P]);
        P = C;
    }
}

template<class Comp>
inline void makeHeap(const Comp& comp, size_t P, device::PKDParticle* particle, size_t N, int dim) {
    if (P >= N)
        return;
    const size_t L = lChild(P);
    const size_t R = rChild(P);
    makeHeap(comp, L, particle, N, dim);
    makeHeap(comp, R, particle, N, dim);
    trickle(comp, P, particle, N, dim);
}
} // namespace megamol::optix_hpg
