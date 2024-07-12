#include "SphereGeometry.h"

#include <sstream>

#include "mmcore/param/BoolParam.h"

#include "sphere.h"

#include "optix/CallGeometry.h"
#include "optix/Utils.h"

#include "optix_stubs.h"


namespace megamol::optix_hpg {
extern "C" const char embedded_sphere_programs[];
//extern "C" const char embedded_sphere_occlusion_programs[];
} // namespace megamol::optix_hpg


megamol::optix_hpg::SphereGeometry::SphereGeometry()
        : _out_geo_slot("outGeo", "")
        , _in_data_slot("inData", "")
        , built_in_intersector_slot_("built-in", "") {
    _out_geo_slot.SetCallback(CallGeometry::ClassName(), CallGeometry::FunctionName(0), &SphereGeometry::get_data_cb);
    _out_geo_slot.SetCallback(
        CallGeometry::ClassName(), CallGeometry::FunctionName(1), &SphereGeometry::get_extents_cb);
    MakeSlotAvailable(&_out_geo_slot);

    _in_data_slot.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&_in_data_slot);

    built_in_intersector_slot_ << new core::param::BoolParam(false);
    MakeSlotAvailable(&built_in_intersector_slot_);
}


megamol::optix_hpg::SphereGeometry::~SphereGeometry() {
    this->Release();
}


bool megamol::optix_hpg::SphereGeometry::create() {
    return true;
}


void megamol::optix_hpg::SphereGeometry::release() {
    if (_geo_buffer != 0) {
        CUDA_CHECK_ERROR(cuMemFree(_geo_buffer));
    }
    /*if (_particle_data != 0) {
        CUDA_CHECK_ERROR(cuMemFree(_particle_data));
    }
    if (color_data_ != 0) {
        CUDA_CHECK_ERROR(cuMemFree(color_data_));
    }*/
    for (auto const& el : particle_data_) {
        CUDA_CHECK_ERROR(cuMemFree(el));
    }
    for (auto const& el : radius_data_) {
        CUDA_CHECK_ERROR(cuMemFree(el));
    }
    for (auto const& el : color_data_) {
        CUDA_CHECK_ERROR(cuMemFree(el));
    }
}


void megamol::optix_hpg::SphereGeometry::init(Context const& ctx) {
    OptixBuiltinISOptions opts = {};
    opts.builtinISModuleType = OPTIX_PRIMITIVE_TYPE_SPHERE;
    opts.usesMotionBlur = false;
    opts.buildFlags = OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
    OPTIX_CHECK_ERROR(optixBuiltinISModuleGet(ctx.GetOptiXContext(), &ctx.GetModuleCompileOptions(),
        &ctx.GetPipelineCompileOptions(), &opts, &sphere_intersector_));

    sphere_module_bi_ = MMOptixModule(embedded_sphere_programs, ctx.GetOptiXContext(), &ctx.GetModuleCompileOptions(),
        &ctx.GetPipelineCompileOptions(), MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
        sphere_intersector_,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "sphere_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "sphere_closesthit"}});
    /*sphere_occlusion_module_bi_ = MMOptixModule(embedded_sphere_occlusion_programs, ctx.GetOptiXContext(),
        &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
        MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP, sphere_intersector_,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "sphere_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "sphere_closesthit_occlusion"}});*/

    sphere_module_ = MMOptixModule(embedded_sphere_programs, ctx.GetOptiXContext(), &ctx.GetModuleCompileOptions(),
        &ctx.GetPipelineCompileOptions(), MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "sphere_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "sphere_closesthit"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "sphere_bounds"}});
    /*sphere_occlusion_module_ = MMOptixModule(embedded_sphere_occlusion_programs, ctx.GetOptiXContext(),
        &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
        MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "sphere_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "sphere_closesthit_occlusion"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "sphere_bounds_occlusion"}});*/

    ++program_version;

    // OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(sphere_module_, &_sbt_record));
}


bool megamol::optix_hpg::SphereGeometry::assertData(geocalls::MultiParticleDataCall& call, Context const& ctx) {
#ifdef MEGAMOL_USE_POWER
    auto power_callbacks = this->frontend_resources.get<frontend_resources::PowerCallbacks>();
#endif

    auto const pl_count = call.GetParticleListCount();

    for (auto const& el : particle_data_) {
        CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
    }
    for (auto const& el : radius_data_) {
        CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
    }
    for (auto const& el : color_data_) {
        CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
    }

    particle_data_.resize(pl_count, 0);
    radius_data_.resize(pl_count, 0);
    color_data_.resize(pl_count, 0);
    std::vector<CUdeviceptr> bounds_data(pl_count);
    std::vector<OptixBuildInput> build_inputs;

#ifdef MEGAMOL_USE_POWER
    size_t total_original_data_size = 0;
#endif

    for (unsigned int pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
        // for now only the first geometry
        auto const& particles = call.AccessParticles(pl_idx);

        auto const p_count = particles.GetCount();
        if (p_count == 0)
            continue;

        /*auto const color_type = particles.GetColourDataType();
        auto const has_color = (color_type != geocalls::SimpleSphericalParticles::COLDATA_NONE) &&
                               (color_type != geocalls::SimpleSphericalParticles::COLDATA_DOUBLE_I) &&
                               (color_type != geocalls::SimpleSphericalParticles::COLDATA_FLOAT_I);
        auto const vert_type = particles.GetVertexDataType();
        auto const has_gobal_radius = vert_type != geocalls::SimpleSphericalParticles::VERTDATA_FLOAT_XYZR;*/

        std::vector<device::Particle> data(p_count);
        auto x_acc = particles.GetParticleStore().GetXAcc();
        auto y_acc = particles.GetParticleStore().GetYAcc();
        auto z_acc = particles.GetParticleStore().GetZAcc();

        auto rad_acc = particles.GetParticleStore().GetRAcc();

        auto cr_acc = particles.GetParticleStore().GetCRAcc();
        auto cg_acc = particles.GetParticleStore().GetCGAcc();
        auto cb_acc = particles.GetParticleStore().GetCBAcc();
        auto ca_acc = particles.GetParticleStore().GetCAAcc();

        for (std::size_t p_idx = 0; p_idx < p_count; ++p_idx) {
            data[p_idx].pos.x = x_acc->Get_f(p_idx);
            data[p_idx].pos.y = y_acc->Get_f(p_idx);
            data[p_idx].pos.z = z_acc->Get_f(p_idx);
        }

        std::vector<float> rad_data;
        if (!has_global_radius(particles)) {
            rad_data.resize(p_count);
            for (std::size_t p_idx = 0; p_idx < p_count; ++p_idx) {
                rad_data[p_idx] = rad_acc->Get_f(p_idx);
            }
        } else {
            rad_data.push_back(particles.GetGlobalRadius());
        }

        auto col_count = p_count;
        if (!has_color(particles)) {
            col_count = 0;
        }
        std::vector<glm::vec4> color_data(col_count);
        if (has_color(particles)) {
            for (std::size_t p_idx = 0; p_idx < col_count; ++p_idx) {
                color_data[p_idx].r = cr_acc->Get_f(p_idx);
                color_data[p_idx].g = cg_acc->Get_f(p_idx);
                color_data[p_idx].b = cb_acc->Get_f(p_idx);
                color_data[p_idx].a = ca_acc->Get_f(p_idx);
            }
            // CUDA_CHECK_ERROR(cuMemFree(color_data_));
            CUDA_CHECK_ERROR(cuMemAllocAsync(&color_data_[pl_idx], col_count * sizeof(glm::vec4), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(
                color_data_[pl_idx], color_data.data(), col_count * sizeof(glm::vec4), ctx.GetExecStream()));
        }
        // CUDA_CHECK_ERROR(cuMemFree(_particle_data));
        CUDA_CHECK_ERROR(
            cuMemAllocAsync(&particle_data_[pl_idx], p_count * sizeof(device::Particle), ctx.GetExecStream()));
        CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(
            particle_data_[pl_idx], data.data(), p_count * sizeof(device::Particle), ctx.GetExecStream()));

        CUDA_CHECK_ERROR(cuMemAllocAsync(&radius_data_[pl_idx], rad_data.size() * sizeof(float), ctx.GetExecStream()));
        CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(
            radius_data_[pl_idx], rad_data.data(), rad_data.size() * sizeof(float), ctx.GetExecStream()));

        unsigned int geo_flag = OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT;
        build_inputs.emplace_back();
        OptixBuildInput& buildInput = build_inputs.back();
        memset(&buildInput, 0, sizeof(OptixBuildInput));

#ifdef MEGAMOL_USE_POWER
        total_original_data_size += data.size() * sizeof(glm::vec3);
#endif

        if (built_in_intersector_slot_.Param<core::param::BoolParam>()->Value()) {
            buildInput.type = OPTIX_BUILD_INPUT_TYPE_SPHERES;
            auto& cp_input = buildInput.sphereArray;
            cp_input.flags = &geo_flag;
            cp_input.numVertices = p_count;
            cp_input.primitiveIndexOffset = 0;
            cp_input.vertexBuffers = &particle_data_[pl_idx];
            cp_input.vertexStrideInBytes = sizeof(device::Particle);
            cp_input.radiusBuffers = &radius_data_[pl_idx];
            cp_input.radiusStrideInBytes = sizeof(float);
            cp_input.numSbtRecords = 1;
            cp_input.sbtIndexOffsetBuffer = NULL;
            cp_input.sbtIndexOffsetSizeInBytes = 0;
            cp_input.sbtIndexOffsetStrideInBytes = 0;
            if (has_global_radius(particles)) {
                cp_input.singleRadius = 1;
            } else {
                cp_input.singleRadius = 0;
            }
        } else {
            CUDA_CHECK_ERROR(
                cuMemAllocAsync(&bounds_data[pl_idx], p_count * sizeof(device::box3f), ctx.GetExecStream()));

            if (has_global_radius(particles)) {
                sphere_module_.ComputeBounds(particle_data_[pl_idx], 0, particles.GetGlobalRadius(),
                    bounds_data[pl_idx], p_count, ctx.GetExecStream());
            } else {
                sphere_module_.ComputeBounds(particle_data_[pl_idx], radius_data_[pl_idx], particles.GetGlobalRadius(),
                    bounds_data[pl_idx], p_count, ctx.GetExecStream());
            }

            //////////////////////////////////////
            // geometry
            //////////////////////////////////////


            buildInput.type = OPTIX_BUILD_INPUT_TYPE_CUSTOM_PRIMITIVES;
            auto& cp_input = buildInput.customPrimitiveArray;
            cp_input.aabbBuffers = &bounds_data[pl_idx];
            cp_input.numPrimitives = p_count;
            cp_input.primitiveIndexOffset = 0;
            cp_input.numSbtRecords = 1;
            cp_input.flags = &geo_flag;
            cp_input.sbtIndexOffsetBuffer = NULL;
            cp_input.sbtIndexOffsetSizeInBytes = 0;
            cp_input.sbtIndexOffsetStrideInBytes = 0;
            cp_input.strideInBytes = 0;
        }

        //SBTRecord<device::SphereGeoData> sbt_record;
        //if (built_in_intersector_slot_.Param<core::param::BoolParam>()->Value()) {
        //    OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(sphere_module_bi_, &sbt_record));
        //} else {
        //    OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(sphere_module_, &sbt_record));
        //}

        //sbt_record.data.particleBufferPtr = (device::Particle*) particle_data_[pl_idx];
        //sbt_record.data.radiusBufferPtr = nullptr;
        //sbt_record.data.colorBufferPtr = nullptr;
        //sbt_record.data.radius = particles.GetGlobalRadius();
        //sbt_record.data.hasGlobalRadius = has_gobal_radius;
        //sbt_record.data.hasColorData = has_color;
        //sbt_record.data.globalColor =
        //    glm::vec4(particles.GetGlobalColour()[0] / 255.f, particles.GetGlobalColour()[1] / 255.f,
        //        particles.GetGlobalColour()[2] / 255.f, particles.GetGlobalColour()[3] / 255.f);

        //if (!has_gobal_radius) {
        //    sbt_record.data.radiusBufferPtr = (float*) radius_data_[pl_idx];
        //}
        //if (has_color) {
        //    sbt_record.data.colorBufferPtr = (glm::vec4*) color_data_[pl_idx];
        //}
        //sbt_records_.push_back(sbt_record);

        //// occlusion stuff
        //SBTRecord<device::SphereGeoData> sbt_record_occlusion;
        //if (built_in_intersector_slot_.Param<core::param::BoolParam>()->Value()) {
        //    OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(sphere_occlusion_module_bi_, &sbt_record_occlusion));
        //} else {
        //    OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(sphere_occlusion_module_, &sbt_record_occlusion));
        //}
        //sbt_record_occlusion.data = sbt_record.data;
        //sbt_records_.push_back(sbt_record_occlusion);

        //++sbt_version;
    }

    OptixAccelBuildOptions accelOptions = {};
    accelOptions.buildFlags = OPTIX_BUILD_FLAG_PREFER_FAST_TRACE | OPTIX_BUILD_FLAG_ALLOW_COMPACTION;
    accelOptions.operation = OPTIX_BUILD_OPERATION_BUILD;
    accelOptions.motionOptions.numKeys = 0;

    OptixAccelBufferSizes bufferSizes = {};
    OPTIX_CHECK_ERROR(optixAccelComputeMemoryUsage(
        ctx.GetOptiXContext(), &accelOptions, build_inputs.data(), build_inputs.size(), &bufferSizes));

    CUdeviceptr geo_temp;
    CUDA_CHECK_ERROR(cuMemFreeAsync(_geo_buffer, ctx.GetExecStream()));
    CUDA_CHECK_ERROR(cuMemAllocAsync(&_geo_buffer, bufferSizes.outputSizeInBytes, ctx.GetExecStream()));
    CUDA_CHECK_ERROR(cuMemAllocAsync(&geo_temp, bufferSizes.tempSizeInBytes, ctx.GetExecStream()));

    CUdeviceptr d_compSize;
    CUDA_CHECK_ERROR(cuMemAllocAsync(&d_compSize, sizeof(std::size_t), ctx.GetExecStream()));
    OptixAccelEmitDesc build_prop = {};
    build_prop.type = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
    build_prop.result = d_compSize;

    OPTIX_CHECK_ERROR(optixAccelBuild(ctx.GetOptiXContext(), ctx.GetExecStream(), &accelOptions, build_inputs.data(),
        build_inputs.size(), geo_temp, bufferSizes.tempSizeInBytes, _geo_buffer, bufferSizes.outputSizeInBytes,
        &_geo_handle, &build_prop, 1));
    CUDA_CHECK_ERROR(cuMemFreeAsync(geo_temp, ctx.GetExecStream()));

    std::size_t compSize;
    CUDA_CHECK_ERROR(cuMemcpyDtoHAsync(&compSize, d_compSize, sizeof(std::size_t), ctx.GetExecStream()));
    CUDA_CHECK_ERROR(cuMemFreeAsync(d_compSize, ctx.GetExecStream()));
#ifdef MEGAMOL_USE_POWER
    if (compSize < bufferSizes.outputSizeInBytes) {
        power_callbacks.add_meta_key_value("GeoSize", std::to_string(compSize));
        core::utility::log::Log::DefaultLog.WriteInfo("[SphereGeometry] Data size BVH only: %d", compSize);
        core::utility::log::Log::DefaultLog.WriteInfo(
            "[SphereGeometry] Temp size BVH only: %d", bufferSizes.tempSizeInBytes);
        core::utility::log::Log::DefaultLog.WriteInfo(
            "[SphereGeometry] Data size with BVH: %d", total_original_data_size + compSize);
    } else {
        power_callbacks.add_meta_key_value("GeoSize", std::to_string(bufferSizes.outputSizeInBytes));
        core::utility::log::Log::DefaultLog.WriteInfo(
            "[SphereGeometry] Data size BVH only: %d", bufferSizes.outputSizeInBytes);
        core::utility::log::Log::DefaultLog.WriteInfo(
            "[SphereGeometry] Temp size BVH only: %d", bufferSizes.tempSizeInBytes);
        core::utility::log::Log::DefaultLog.WriteInfo(
            "[SphereGeometry] Data size with BVH: %d", total_original_data_size + bufferSizes.outputSizeInBytes);
    }
    power_callbacks.add_meta_key_value("OriginalDataSize", std::to_string(total_original_data_size));
#endif
    if (compSize < bufferSizes.outputSizeInBytes) {
        CUdeviceptr comp_geo_buffer;
        CUDA_CHECK_ERROR(cuMemAllocAsync(&comp_geo_buffer, compSize, ctx.GetExecStream()));
        OptixTraversableHandle compacted_geo_handle = 0;
        OPTIX_CHECK_ERROR(optixAccelCompact(
            ctx.GetOptiXContext(), ctx.GetExecStream(), _geo_handle, comp_geo_buffer, compSize, &compacted_geo_handle));

        auto temp_handle = _geo_handle;
        auto temp_buffer = _geo_buffer;
        _geo_handle = compacted_geo_handle;
        _geo_buffer = comp_geo_buffer;

        CUDA_CHECK_ERROR(cuMemFreeAsync(temp_buffer, ctx.GetExecStream()));
    }

    ++geo_version;

    // CUDA_CHECK_ERROR(cuMemFree(bounds_data));
    if (!built_in_intersector_slot_.Param<core::param::BoolParam>()->Value()) {
        for (auto const& el : bounds_data) {
            CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
        }
    }

    //////////////////////////////////////
    // end geometry
    //////////////////////////////////////

    return true;
}

bool megamol::optix_hpg::SphereGeometry::createSBTRecords(geocalls::MultiParticleDataCall& call, Context const& ctx) {
    auto const pl_count = call.GetParticleListCount();

    sbt_records_.clear();

    for (unsigned int pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
        auto const& particles = call.AccessParticles(pl_idx);

        auto const p_count = particles.GetCount();
        if (p_count == 0)
            continue;

        SBTRecord<device::SphereGeoData> sbt_record;
        if (built_in_intersector_slot_.Param<core::param::BoolParam>()->Value()) {
            OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(sphere_module_bi_, &sbt_record));
        } else {
            OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(sphere_module_, &sbt_record));
        }

        sbt_record.data.particleBufferPtr = (device::Particle*) particle_data_[pl_idx];
        sbt_record.data.radiusBufferPtr = nullptr;
        sbt_record.data.colorBufferPtr = nullptr;
        sbt_record.data.radius = particles.GetGlobalRadius();
        sbt_record.data.hasGlobalRadius = has_global_radius(particles);
        sbt_record.data.hasColorData = has_color(particles);
        sbt_record.data.globalColor =
            glm::vec4(particles.GetGlobalColour()[0] / 255.f, particles.GetGlobalColour()[1] / 255.f,
                particles.GetGlobalColour()[2] / 255.f, particles.GetGlobalColour()[3] / 255.f);

        if (!has_global_radius(particles)) {
            sbt_record.data.radiusBufferPtr = (float*) radius_data_[pl_idx];
        }
        if (has_color(particles)) {
            sbt_record.data.colorBufferPtr = (glm::vec4*) color_data_[pl_idx];
        }
        sbt_records_.push_back(sbt_record);

        // occlusion stuff
        /*SBTRecord<device::SphereGeoData> sbt_record_occlusion;
        if (built_in_intersector_slot_.Param<core::param::BoolParam>()->Value()) {
            OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(sphere_occlusion_module_bi_, &sbt_record_occlusion));
        } else {
            OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(sphere_occlusion_module_, &sbt_record_occlusion));
        }
        sbt_record_occlusion.data = sbt_record.data;
        sbt_records_.push_back(sbt_record_occlusion);*/
    }

    ++sbt_version;

    return true;
}


bool megamol::optix_hpg::SphereGeometry::get_data_cb(core::Call& c) {
    auto out_geo = dynamic_cast<CallGeometry*>(&c);
    if (out_geo == nullptr)
        return false;
    auto in_data = _in_data_slot.CallAs<geocalls::MultiParticleDataCall>();
    if (in_data == nullptr)
        return false;

    auto const ctx = out_geo->get_ctx();

    static bool not_init = true;
    if (not_init) {
        init(*ctx);
        not_init = false;
    }

    in_data->SetFrameID(out_geo->FrameID());
    if (!(*in_data)(1))
        return false;
    if (!(*in_data)(0))
        return false;

    if (in_data->FrameID() != _frame_id || in_data->DataHash() != _data_hash || built_in_intersector_slot_.IsDirty()) {
        if (!assertData(*in_data, *ctx))
            return false;
        createSBTRecords(*in_data, *ctx);
        if (built_in_intersector_slot_.IsDirty()) {
            ++program_version;
            built_in_intersector_slot_.ResetDirty();
        }
        _frame_id = in_data->FrameID();
        _data_hash = in_data->DataHash();
        
    }
    
    if (built_in_intersector_slot_.Param<core::param::BoolParam>()->Value()) {
        program_groups_[0] = sphere_module_bi_;
        //program_groups_[1] = sphere_occlusion_module_bi_;
    } else {
        program_groups_[0] = sphere_module_;
        //program_groups_[1] = sphere_occlusion_module_;
    }

    out_geo->set_handle(&_geo_handle, geo_version);
    out_geo->set_program_groups(program_groups_.data(), program_groups_.size(), program_version);
    out_geo->set_record(
        sbt_records_.data(), sbt_records_.size(), sizeof(SBTRecord<device::SphereGeoData>), sbt_version);

    return true;
}


bool megamol::optix_hpg::SphereGeometry::get_extents_cb(core::Call& c) {
    auto out_geo = dynamic_cast<CallGeometry*>(&c);
    if (out_geo == nullptr)
        return false;
    auto in_data = _in_data_slot.CallAs<geocalls::MultiParticleDataCall>();
    if (in_data == nullptr)
        return false;

    if ((*in_data)(1)) {
        out_geo->SetFrameCount(in_data->FrameCount());
        out_geo->AccessBoundingBoxes() = in_data->AccessBoundingBoxes();
    }

    return true;
}
