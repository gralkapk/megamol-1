#include "PKDGeometry.h"

#include "mmcore/param/BoolParam.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/IntParam.h"

#include <glm/glm.hpp>

#include <omp.h>

#include <tbb/parallel_for.h>

#include "CallGeometry.h"

namespace megamol::optix_hpg {
extern "C" const char embedded_pkd_programs[];

int arg_max(glm::vec3 const& v) {
    int biggestDim = 0;
    for (int i = 1; i < 3; i++)
        if ((v[i]) > (v[biggestDim]))
            biggestDim = i;
    return biggestDim;
}

void recBuild(size_t /* root node */ P, device::PKDParticle* particle, size_t N, box3f bounds) {
    if (P >= N)
        return;

    int dim = arg_max(bounds.upper - bounds.lower);

    const size_t L = lChild(P);
    const size_t R = rChild(P);
    const bool lValid = (L < N);
    const bool rValid = (R < N);
    makeHeap(std::greater<float>(), L, particle, N, dim);
    makeHeap(std::less<float>(), R, particle, N, dim);

    if (rValid) {
        while (particle[L].pos[dim] > particle[R].pos[dim]) {
            std::swap(particle[L], particle[R]);
            trickle(std::greater<float>(), L, particle, N, dim);
            trickle(std::less<float>(), R, particle, N, dim);
        }
        if (particle[L].pos[dim] > particle[P].pos[dim]) {
            std::swap(particle[L], particle[P]);
            particle[L].dim = dim;
        } else if (particle[R].pos[dim] < particle[P].pos[dim]) {
            std::swap(particle[R], particle[P]);
            particle[R].dim = dim;
        } else
            /* nothing, root fits */;
    } else if (lValid) {
        if (particle[L].pos[dim] > particle[P].pos[dim]) {
            std::swap(particle[L], particle[P]);
            particle[L].dim = dim;
        }
    }

    box3f lBounds = bounds;
    box3f rBounds = bounds;
    lBounds.upper[dim] = rBounds.lower[dim] = particle[P].pos[dim];
    particle[P].dim = dim;

    tbb::parallel_for(0, 2, [&](int childID) {
        if (childID) {
            recBuild(L, particle, N, lBounds);
        } else {
            recBuild(R, particle, N, rBounds);
        }
    });
}

void makePKD(std::vector<device::PKDParticle>& particles, box3f bounds) {
    recBuild(/*node:*/ 0, particles.data(), particles.size(), bounds);
}

PKDGeometry::PKDGeometry()
        : out_geo_slot_("outGeo", "")
        , in_data_slot_("inData", "")
        , mode_slot_("mode", "")
        , compression_slot_("compression", "")
        , threshold_slot_("threshold", "") {
    out_geo_slot_.SetCallback(CallGeometry::ClassName(), CallGeometry::FunctionName(0), &PKDGeometry::get_data_cb);
    out_geo_slot_.SetCallback(CallGeometry::ClassName(), CallGeometry::FunctionName(1), &PKDGeometry::get_extents_cb);
    MakeSlotAvailable(&out_geo_slot_);

    in_data_slot_.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&in_data_slot_);

    auto ep = new core::param::EnumParam(static_cast<int>(PKDMode::STANDARD));
    ep->SetTypePair(static_cast<int>(PKDMode::STANDARD), "Standard");
    ep->SetTypePair(static_cast<int>(PKDMode::TREELETS), "Treelets");
    mode_slot_ << ep;
    MakeSlotAvailable(&mode_slot_);

    compression_slot_ << new core::param::BoolParam(false);
    MakeSlotAvailable(&compression_slot_);

    threshold_slot_ << new core::param::IntParam(256, 16);
    MakeSlotAvailable(&threshold_slot_);
}

PKDGeometry::~PKDGeometry() {
    this->Release();
}

bool PKDGeometry::create() {
    return true;
}

void PKDGeometry::release() {
    for (auto& el : particle_data_) {
        CUDA_CHECK_ERROR(cuMemFree(el));
    }
    for (auto& el : radius_data_) {
        CUDA_CHECK_ERROR(cuMemFree(el));
    }
    for (auto& el : color_data_) {
        CUDA_CHECK_ERROR(cuMemFree(el));
    }
}

bool PKDGeometry::get_data_cb(core::Call& c) {
    auto out_geo = dynamic_cast<CallGeometry*>(&c);
    if (out_geo == nullptr)
        return false;
    auto in_data = in_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
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

    if (in_data->FrameID() != frame_id_ || in_data->DataHash() != data_hash_) {
        if (!assert_data(*in_data, *ctx))
            return false;
        createSBTRecords(*in_data, *ctx);
        frame_id_ = in_data->FrameID();
        data_hash_ = in_data->DataHash();
    }
    program_groups_[0] = pkd_module_;
    program_groups_[1] = pkd_module_occlusion_;
    

    out_geo->set_handle(&geo_handle_);
    out_geo->set_program_groups(program_groups_.data(), program_groups_.size(), program_version);
    out_geo->set_record(
        sbt_records_.data(), sbt_records_.size(), sizeof(SBTRecord<device::PKDGeoData>), sbt_version);

    return true;
}

bool PKDGeometry::get_extents_cb(core::Call& c) {
    auto out_geo = dynamic_cast<CallGeometry*>(&c);
    if (out_geo == nullptr)
        return false;
    auto in_data = in_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (in_data == nullptr)
        return false;

    if ((*in_data)(1)) {
        out_geo->SetFrameCount(in_data->FrameCount());
        out_geo->AccessBoundingBoxes() = in_data->AccessBoundingBoxes();
    }

    return true;
}

bool PKDGeometry::init(Context const& ctx) {
    pkd_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(), &ctx.GetModuleCompileOptions(),
        &ctx.GetPipelineCompileOptions(), MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "pkd_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "pkd_closesthit"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "pkd_bounds"}});

    pkd_module_occlusion_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(), &ctx.GetModuleCompileOptions(),
        &ctx.GetPipelineCompileOptions(), MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "pkd_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "pkd_closesthit_occlusion"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "pkd_bounds"}});

    ++program_version;

    return true;
}

bool PKDGeometry::assert_data(geocalls::MultiParticleDataCall const& call, Context const& ctx) {
    free_data(ctx);

    auto const pl_count = call.GetParticleListCount();

    particle_data_.resize(pl_count, 0);
    radius_data_.resize(pl_count, 0);
    color_data_.resize(pl_count, 0);
    local_boxes_.resize(
        pl_count, std::make_pair(glm::vec3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                                     std::numeric_limits<float>::max()),
                      glm::vec3(std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(),
                          std::numeric_limits<float>::lowest())));
    std::vector<CUdeviceptr> bounds_data(pl_count);
    std::vector<OptixBuildInput> build_inputs;

    for (unsigned int pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
        auto const& particles = call.AccessParticles(pl_idx);

        auto const p_count = particles.GetCount();
        if (p_count == 0)
            continue;

        std::vector<device::PKDParticle> data(p_count);
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

        // TODO Compute tight bounds
        auto const max_threads = omp_get_max_threads();
        glm::vec3 initial_lower = glm::vec3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                      std::numeric_limits<float>::max()),
                  initial_upper = glm::vec3(std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(),
                      std::numeric_limits<float>::lowest());
        std::vector<std::pair<glm::vec3, glm::vec3>> local_boxes(
            max_threads, std::make_pair(initial_lower, initial_upper));
#pragma omp parallel for shared(local_boxes)
        for (int64_t p_idx = 0; p_idx < p_count; ++p_idx) {
            auto const thread_num = omp_get_thread_num();
            auto& box = local_boxes[thread_num];
            glm::vec3 pos(x_acc->Get_f(p_idx), y_acc->Get_f(p_idx), z_acc->Get_f(p_idx));
            auto const rad = rad_acc->Get_f(p_idx);
            auto const new_lower = pos - rad;
            auto const new_upper = pos + rad;
            std::get<0>(box).x = std::min(std::get<0>(box).x, new_lower.x);
            std::get<0>(box).y = std::min(std::get<0>(box).y, new_lower.y);
            std::get<0>(box).z = std::min(std::get<0>(box).z, new_lower.z);
            std::get<1>(box).x = std::max(std::get<1>(box).x, new_upper.x);
            std::get<1>(box).y = std::max(std::get<1>(box).y, new_upper.y);
            std::get<1>(box).z = std::max(std::get<1>(box).z, new_upper.z);
        }

        std::array<float, 6> local_box = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
            std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest(),
            std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};
        for (auto const& el : local_boxes) {
            local_box[0] = std::min(std::get<0>(el).x, local_box[0]);
            local_box[1] = std::min(std::get<0>(el).y, local_box[1]);
            local_box[2] = std::min(std::get<0>(el).z, local_box[2]);
            local_box[3] = std::max(std::get<1>(el).x, local_box[3]);
            local_box[4] = std::max(std::get<1>(el).y, local_box[4]);
            local_box[5] = std::max(std::get<1>(el).z, local_box[5]);
        }

        box3f op_box;
        op_box.lower = glm::vec3(local_box[0], local_box[1], local_box[2]);
        op_box.upper = glm::vec3(local_box[3], local_box[4], local_box[5]);
        makePKD(data, op_box);
        local_boxes_[pl_idx] = std::make_pair(op_box.lower, op_box.upper);


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
            CUDA_CHECK_ERROR(cuMemAllocAsync(&color_data_[pl_idx], col_count * sizeof(glm::vec4), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(
                color_data_[pl_idx], color_data.data(), col_count * sizeof(glm::vec4), ctx.GetExecStream()));
        }
        CUDA_CHECK_ERROR(
            cuMemAllocAsync(&particle_data_[pl_idx], p_count * sizeof(device::PKDParticle), ctx.GetExecStream()));
        CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(
            particle_data_[pl_idx], data.data(), p_count * sizeof(device::PKDParticle), ctx.GetExecStream()));

        CUDA_CHECK_ERROR(cuMemAllocAsync(&radius_data_[pl_idx], rad_data.size() * sizeof(float), ctx.GetExecStream()));
        CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(
            radius_data_[pl_idx], rad_data.data(), rad_data.size() * sizeof(float), ctx.GetExecStream()));

        unsigned int geo_flag = OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT;
        build_inputs.emplace_back();
        OptixBuildInput& buildInput = build_inputs.back();
        memset(&buildInput, 0, sizeof(OptixBuildInput));


        CUDA_CHECK_ERROR(cuMemAllocAsync(&bounds_data[pl_idx], 1 * sizeof(box3f), ctx.GetExecStream()));
        CUDA_CHECK_ERROR(
            cuMemcpyHtoDAsync(bounds_data[pl_idx], local_box.data(), 6 * sizeof(float), ctx.GetExecStream()));


        //////////////////////////////////////
        // geometry
        //////////////////////////////////////


        buildInput.type = OPTIX_BUILD_INPUT_TYPE_CUSTOM_PRIMITIVES;
        auto& cp_input = buildInput.customPrimitiveArray;
        cp_input.aabbBuffers = &bounds_data[pl_idx];
        cp_input.numPrimitives = 1;
        cp_input.primitiveIndexOffset = 0;
        cp_input.numSbtRecords = 1;
        cp_input.flags = &geo_flag;
        cp_input.sbtIndexOffsetBuffer = NULL;
        cp_input.sbtIndexOffsetSizeInBytes = 0;
        cp_input.sbtIndexOffsetStrideInBytes = 0;
        cp_input.strideInBytes = 0;
    }

    OptixAccelBuildOptions accelOptions = {};
    accelOptions.buildFlags = OPTIX_BUILD_FLAG_PREFER_FAST_TRACE | OPTIX_BUILD_FLAG_ALLOW_COMPACTION;
    accelOptions.operation = OPTIX_BUILD_OPERATION_BUILD;
    accelOptions.motionOptions.numKeys = 0;

    OptixAccelBufferSizes bufferSizes = {};
    OPTIX_CHECK_ERROR(optixAccelComputeMemoryUsage(
        ctx.GetOptiXContext(), &accelOptions, build_inputs.data(), build_inputs.size(), &bufferSizes));

    CUdeviceptr geo_temp;
    CUDA_CHECK_ERROR(cuMemFreeAsync(geo_buffer_, ctx.GetExecStream()));
    CUDA_CHECK_ERROR(cuMemAllocAsync(&geo_buffer_, bufferSizes.outputSizeInBytes, ctx.GetExecStream()));
    CUDA_CHECK_ERROR(cuMemAllocAsync(&geo_temp, bufferSizes.tempSizeInBytes, ctx.GetExecStream()));

    CUdeviceptr d_compSize;
    CUDA_CHECK_ERROR(cuMemAllocAsync(&d_compSize, sizeof(std::size_t), ctx.GetExecStream()));
    OptixAccelEmitDesc build_prop = {};
    build_prop.type = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
    build_prop.result = d_compSize;

    OPTIX_CHECK_ERROR(optixAccelBuild(ctx.GetOptiXContext(), ctx.GetExecStream(), &accelOptions, build_inputs.data(),
        build_inputs.size(), geo_temp, bufferSizes.tempSizeInBytes, geo_buffer_, bufferSizes.outputSizeInBytes,
        &geo_handle_, &build_prop, 1));

    std::size_t compSize;
    CUDA_CHECK_ERROR(cuMemcpyDtoHAsync(&compSize, d_compSize, sizeof(std::size_t), ctx.GetExecStream()));
    CUDA_CHECK_ERROR(cuMemFreeAsync(d_compSize, ctx.GetExecStream()));
    if (compSize < bufferSizes.outputSizeInBytes) {
        CUdeviceptr comp_geo_buffer;
        CUDA_CHECK_ERROR(cuMemAllocAsync(&comp_geo_buffer, compSize, ctx.GetExecStream()));
        OptixTraversableHandle compacted_geo_handle = 0;
        OPTIX_CHECK_ERROR(optixAccelCompact(
            ctx.GetOptiXContext(), ctx.GetExecStream(), geo_handle_, comp_geo_buffer, compSize, &compacted_geo_handle));

        auto temp_handle = geo_handle_;
        auto temp_buffer = geo_buffer_;
        geo_handle_ = compacted_geo_handle;
        geo_buffer_ = comp_geo_buffer;

        CUDA_CHECK_ERROR(cuMemFreeAsync(temp_buffer, ctx.GetExecStream()));
    }

    CUDA_CHECK_ERROR(cuMemFreeAsync(geo_temp, ctx.GetExecStream()));
    for (auto const& el : bounds_data) {
        CUDA_CHECK_ERROR(cuMemFreeAsync(el, ctx.GetExecStream()));
    }

    //////////////////////////////////////
    // end geometry
    //////////////////////////////////////

    return true;
}

bool PKDGeometry::createSBTRecords(geocalls::MultiParticleDataCall const& call, Context const& ctx) {
    auto const pl_count = call.GetParticleListCount();

    sbt_records_.clear();

    for (unsigned int pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
        auto const& particles = call.AccessParticles(pl_idx);

        auto const p_count = particles.GetCount();
        if (p_count == 0)
            continue;

        SBTRecord<device::PKDGeoData> sbt_record;
        OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(pkd_module_, &sbt_record));

        sbt_record.data.particleBufferPtr = (device::PKDParticle*) particle_data_[pl_idx];
        sbt_record.data.radiusBufferPtr = nullptr;
        sbt_record.data.colorBufferPtr = nullptr;
        sbt_record.data.radius = particles.GetGlobalRadius();
        sbt_record.data.hasGlobalRadius = has_global_radius(particles);
        sbt_record.data.hasColorData = has_color(particles);
        sbt_record.data.globalColor =
            glm::vec4(particles.GetGlobalColour()[0] / 255.f, particles.GetGlobalColour()[1] / 255.f,
                particles.GetGlobalColour()[2] / 255.f, particles.GetGlobalColour()[3] / 255.f);
        sbt_record.data.particleCount = p_count;
        sbt_record.data.worldBounds_lower = std::get<0>(local_boxes_[pl_idx]);
        sbt_record.data.worldBounds_upper = std::get<1>(local_boxes_[pl_idx]);

        if (!has_global_radius(particles)) {
            sbt_record.data.radiusBufferPtr = (float*) radius_data_[pl_idx];
        }
        if (has_color(particles)) {
            sbt_record.data.colorBufferPtr = (glm::vec4*) color_data_[pl_idx];
        }
        sbt_records_.push_back(sbt_record);

        // occlusion stuff
        SBTRecord<device::PKDGeoData> sbt_record_occlusion;
        OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(pkd_module_occlusion_, &sbt_record_occlusion));

        sbt_record_occlusion.data = sbt_record.data;
        sbt_records_.push_back(sbt_record_occlusion);
    }

    ++sbt_version;

    return true;
}

} // namespace megamol::optix_hpg
