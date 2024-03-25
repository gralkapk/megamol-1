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

#include "TreeletCache.h"

#include "nvcomp.hpp"
#include "nvcomp/nvcompManagerFactory.hpp"

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
        for (auto& pel : comp_particle_data_) {
            for (auto& el : pel) {
                CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
            }
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

    void process_treelet_requests(int pl_idx, std::vector<int> const& treelet_reqs);

    std::shared_ptr<TreeletCache> treelet_cache_;
    std::shared_ptr<nvcomp::GdeflateManager> coder_;

    core::CalleeSlot out_geo_slot_;

    core::CallerSlot in_data_slot_;

    core::param::ParamSlot mode_slot_;

    core::param::ParamSlot compression_slot_;

    core::param::ParamSlot threshold_slot_;

    core::param::ParamSlot entropy_slot_;

    std::vector<SBTRecord<device::PKDGeoData>> sbt_records_;

    std::vector<SBTRecord<device::TreeletsGeoData>> treelets_sbt_records_;

    std::vector<SBTRecord<device::QTreeletsGeoData>> comp_treelets_sbt_records_;

    std::array<OptixProgramGroup, 2> program_groups_;

    std::vector<CUdeviceptr> particle_data_;

    std::vector<std::vector<CUdeviceptr>> comp_particle_data_;

    //std::vector<CUdeviceptr> radius_data_;

    std::vector<CUdeviceptr> color_data_;

    std::vector<CUdeviceptr> treelets_data_;

    std::vector<device::box3f> local_boxes_;

    CUdeviceptr geo_buffer_ = 0;

    OptixTraversableHandle geo_handle_;

    MMOptixModule pkd_module_;

    MMOptixModule pkd_occlusion_module_;

    MMOptixModule treelets_module_;

    MMOptixModule treelets_occlusion_module_;

    MMOptixModule comp_treelets_module_;

    MMOptixModule comp_treelets_occlusion_module_;

    uint64_t sbt_version = 0;

    uint64_t program_version = 0;

    uint64_t geo_version = 0;

    unsigned int frame_id_ = std::numeric_limits<unsigned int>::max();

    std::size_t data_hash_ = std::numeric_limits<std::size_t>::max();
};
} // namespace megamol::optix_hpg
