#pragma once

#include "geometry_calls/MultiParticleDataCall.h"
#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"

#include "mmcore/param/ParamSlot.h"

#include "optix/CallContext.h"

#include "cuda.h"

#include "SBT.h"

#include "sphere.h"

#include "MMOptixModule.h"

#include "optix/Context.h"

#ifdef MEGAMOL_USE_POWER
#include "PowerCallbacks.h"
#endif

namespace megamol::optix_hpg {
class SphereGeometry : public core::Module {
public:
    static void requested_lifetime_resources(frontend_resources::ResourceRequest& req) {
        core::Module::requested_lifetime_resources(req);
#ifdef MEGAMOL_USE_POWER
        req.require<frontend_resources::PowerCallbacks>();
#endif
    }

    static const char* ClassName(void) {
        return "SphereGeometry";
    }

    static const char* Description(void) {
        return "Sphere Geometry for OptiX";
    }

    static bool IsAvailable(void) {
        return true;
    }

    SphereGeometry();

    virtual ~SphereGeometry();

protected:
    bool create() override;

    void release() override;

private:
    void init(Context const& ctx);

    bool assertData(geocalls::MultiParticleDataCall& call, Context const& ctx);

    bool createSBTRecords(geocalls::MultiParticleDataCall& call, Context const& ctx);

    bool get_data_cb(core::Call& c);

    bool get_extents_cb(core::Call& c);

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

    core::CalleeSlot _out_geo_slot;

    core::CallerSlot _in_data_slot;

#if OPTIX_VERSION >= 80000
    core::param::ParamSlot built_in_intersector_slot_;
    MMOptixModule sphere_module_bi_;
#endif

    MMOptixModule sphere_module_;

    //MMOptixModule sphere_occlusion_module_;

    //MMOptixModule sphere_occlusion_module_bi_;

    OptixModule sphere_intersector_;

    std::vector<SBTRecord<device::SphereGeoData>> sbt_records_;

    std::array<OptixProgramGroup, 1> program_groups_;

    CUdeviceptr _geo_buffer = 0;

    std::vector<CUdeviceptr> particle_data_;

    std::vector<CUdeviceptr> radius_data_;

    std::vector<CUdeviceptr> color_data_;

    OptixTraversableHandle _geo_handle;

    unsigned int _frame_id = std::numeric_limits<unsigned int>::max();

    std::size_t _data_hash = std::numeric_limits<std::size_t>::max();

    uint64_t sbt_version = 0;

    uint64_t program_version = 0;

    uint64_t geo_version = 0;
};
} // namespace megamol::optix_hpg
