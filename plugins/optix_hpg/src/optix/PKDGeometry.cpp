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

// BEGIN PKD

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

void makePKD(std::vector<device::PKDParticle>& particles, size_t begin, size_t end, box3f bounds) {
    recBuild(/*node:*/ 0, particles.data() + begin, end - begin, bounds);
}

// END PKD

// BEGIN TREELETS

inline size_t sort_partition(
    std::vector<device::PKDParticle>& particles, size_t begin, size_t end, box3f bounds, int& splitDim) {
    // -------------------------------------------------------
    // determine split pos
    // -------------------------------------------------------
    splitDim = arg_max(bounds.upper - bounds.lower);
    //float splitPos = bounds.center()[splitDim];
    float splitPos = (0.5f * (bounds.upper + bounds.lower))[splitDim];

    // -------------------------------------------------------
    // now partition ...
    // -------------------------------------------------------
    size_t mid = begin;
    size_t l = begin, r = (end - 1);
    // quicksort partition:
    while (l <= r) {
        while (l < r && particles[l].pos[splitDim] < splitPos)
            ++l;
        while (l < r && particles[r].pos[splitDim] >= splitPos)
            --r;
        if (l == r) {
            mid = l;
            break;
        }

        std::swap(particles[l], particles[r]);
    }

    // catch-all for extreme cases where all particles are on the same
    // spot, and can't be split:
    if (mid == begin || mid == end)
        mid = (begin + end) / 2;

    return mid;
}

box3f extendBounds(
    std::vector<device::PKDParticle> const& particles, size_t begin, size_t end, float radius) {
    box3f bounds;
    for (int64_t p_idx = begin; p_idx < end; ++p_idx) {
        auto const new_lower = particles[p_idx].pos - radius;
        auto const new_upper = particles[p_idx].pos + radius;
        bounds.extend(new_lower);
        bounds.extend(new_upper);
    }

    return bounds;
}

template<typename MakeLeafLambda>
void partitionRecursively(
    std::vector<device::PKDParticle>& particles, size_t begin, size_t end, const MakeLeafLambda& makeLeaf) {
    if (makeLeaf(begin, end, false))
        // could make into a leaf, done.
        return;

    // -------------------------------------------------------
    // parallel bounding box computation
    // -------------------------------------------------------
    box3f bounds;
    std::mutex boundsMutex;

    const size_t blockSize = 32 * 1024;
    const size_t numTasks = end - begin;
    const size_t numBlocks = (numTasks + blockSize - 1) / blockSize;
    /*parallel_for(numBlocks, [&](size_t blockID) {
        size_t block_begin = begin + blockID * blockSize;
        taskFunction(block_begin, std::min(block_begin + blockSize, end));
    });*/

    auto task = [&](size_t blockBegin, size_t blockEnd) {
        //box3f blockBounds;
        /*for (size_t i = blockBegin; i < blockEnd; i++)
            blockBounds.extend(model->particles[i].pos);*/
        auto const blockBounds = extendBounds(particles, blockBegin, blockEnd, 0.f);
        std::lock_guard<std::mutex> lock(boundsMutex);
        bounds.extend(blockBounds);
        /*bounds.lower.x = std::min(bounds.lower.x, blockBounds.lower.x);
        bounds.lower.y = std::min(bounds.lower.y, blockBounds.lower.y);
        bounds.lower.z = std::min(bounds.lower.z, blockBounds.lower.z);
        bounds.upper.x = std::max(bounds.upper.x, blockBounds.upper.x);
        bounds.upper.y = std::max(bounds.upper.y, blockBounds.upper.y);
        bounds.upper.z = std::max(bounds.upper.z, blockBounds.upper.z);*/
    };
    tbb::parallel_for((size_t) 0, numBlocks, [&](size_t blockID) {
        size_t block_begin = begin + blockID * blockSize;
        task(block_begin, std::min(block_begin + blockSize, end));
    });

    //parallel_for_blocked(begin, end, 32 * 1024, [&](size_t blockBegin, size_t blockEnd) {
    //    box3f blockBounds;
    //    /*for (size_t i = blockBegin; i < blockEnd; i++)
    //        blockBounds.extend(model->particles[i].pos);*/
    //    std::tie(blockBounds.lower, blockBounds.upper) = extendBounds(particles, blockBegin, blockEnd, 0.f);
    //    std::lock_guard<std::mutex> lock(boundsMutex);
    //    /*bounds.extend(blockBounds);*/
    //    bounds.lower.x = std::min(bounds.lower.x, blockBounds.lower.x);
    //    bounds.lower.y = std::min(bounds.lower.y, blockBounds.lower.y);
    //    bounds.lower.z = std::min(bounds.lower.z, blockBounds.lower.z);
    //    bounds.upper.x = std::max(bounds.upper.x, blockBounds.upper.x);
    //    bounds.upper.y = std::max(bounds.upper.y, blockBounds.upper.y);
    //    bounds.upper.z = std::max(bounds.upper.z, blockBounds.upper.z);
    //});

    int splitDim;
    auto mid = sort_partition(particles, begin, end, bounds, splitDim);

    // -------------------------------------------------------
    // and recurse ...
    // -------------------------------------------------------
    tbb::parallel_for(0, 2, [&](int side) {
        if (side)
            partitionRecursively(particles, begin, mid, makeLeaf);
        else
            partitionRecursively(particles, mid, end, makeLeaf);
    });
}

std::vector<device::PKDlet> prePartition_inPlace(
    std::vector<device::PKDParticle>& particles, size_t maxSize, float radius) {
    std::mutex resultMutex;
    std::vector<device::PKDlet> result;

    partitionRecursively(particles, 0ULL, particles.size(), [&](size_t begin, size_t end, bool force) {
        /*bool makeLeaf() :*/
        const size_t size = end - begin;
        if (size > maxSize && !force)
            return false;

        device::PKDlet treelet;
        treelet.begin = begin;
        treelet.end = end;
        /*treelet.bounds = box3f();
        for (size_t i = begin; i < end; ++i) {
            treelet.bounds.extend(model->particles[i].pos - radius);
            treelet.bounds.extend(model->particles[i].pos + radius);
        }*/
        treelet.bounds = extendBounds(particles, begin, end, radius);

        std::lock_guard<std::mutex> lock(resultMutex);
        result.push_back(treelet);
        return true;
    });

    return std::move(result);
}

// END TREELETS

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
    for (auto& el : treelets_data_) {
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

    if (in_data->FrameID() != frame_id_ || in_data->DataHash() != data_hash_ || mode_slot_.IsDirty()) {
        if (!assert_data(*in_data, *ctx))
            return false;
        createSBTRecords(*in_data, *ctx);
        frame_id_ = in_data->FrameID();
        data_hash_ = in_data->DataHash();
        if (mode_slot_.IsDirty()) {
            ++program_version;
        }
        mode_slot_.ResetDirty();
    }
    
    if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) {
        program_groups_[0] = treelets_module_;
        program_groups_[1] = treelets_occlusion_module_;
    } else {
        program_groups_[0] = pkd_module_;
        program_groups_[1] = pkd_module_occlusion_;
    }


    out_geo->set_handle(&geo_handle_, geo_version);
    out_geo->set_program_groups(program_groups_.data(), program_groups_.size(), program_version);
    if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) {
        out_geo->set_record(treelets_sbt_records_.data(), treelets_sbt_records_.size(),
            sizeof(SBTRecord<device::TreeletsGeoData>), sbt_version);
    } else {
        out_geo->set_record(
            sbt_records_.data(), sbt_records_.size(), sizeof(SBTRecord<device::PKDGeoData>), sbt_version);
    }

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
    // TODO Modules for treelets

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

    treelets_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(), &ctx.GetModuleCompileOptions(),
        &ctx.GetPipelineCompileOptions(), MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelets_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "treelets_closesthit"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});

    treelets_occlusion_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
        &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
        MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelets_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "treelets_closesthit_occlusion"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});

    ++program_version;

    return true;
}

bool PKDGeometry::assert_data(geocalls::MultiParticleDataCall const& call, Context const& ctx) {
    free_data(ctx);

    auto const pl_count = call.GetParticleListCount();

    particle_data_.resize(pl_count, 0);
    radius_data_.resize(pl_count, 0);
    color_data_.resize(pl_count, 0);
    treelets_data_.resize(pl_count, 0);
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
        if (p_count == 0 || !has_global_radius(particles)) {
            if (!has_global_radius(particles)) {
                core::utility::log::Log::DefaultLog.WriteWarn("[PKDGeometry]: Per-particle radius not supported");
            }
            continue;
        }

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

        std::array<float, 6> local_box = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
            std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest(),
            std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};

        std::vector<device::PKDlet> treelets;
        if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) {
            // TODO separate particles into a set of treelets
            treelets = prePartition_inPlace(
                data, threshold_slot_.Param<core::param::IntParam>()->Value(), particles.GetGlobalRadius());

            tbb::parallel_for((size_t) 0, treelets.size(), [&](size_t treeletID) {
                /*box3f bounds;
                bounds.lower = treelets[treeletID].bounds_lower;
                bounds.upper = treelets[treeletID].bounds_upper;*/
                makePKD(data, treelets[treeletID].begin, treelets[treeletID].end, treelets[treeletID].bounds);
            });

            CUDA_CHECK_ERROR(cuMemAllocAsync(
                &treelets_data_[pl_idx], treelets.size() * sizeof(device::PKDlet), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(treelets_data_[pl_idx], treelets.data(),
                treelets.size() * sizeof(device::PKDlet), ctx.GetExecStream()));

            // TODO compress data if requested
        } else {
            auto const max_threads = omp_get_max_threads();
            glm::vec3 initial_lower = glm::vec3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                          std::numeric_limits<float>::max()),
                      initial_upper = glm::vec3(std::numeric_limits<float>::lowest(),
                          std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest());
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

        if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) {
            // TODO set of treelet boxes
            CUDA_CHECK_ERROR(
                cuMemAllocAsync(&bounds_data[pl_idx], treelets.size() * sizeof(box3f), ctx.GetExecStream()));
            std::vector<box3f> treelet_boxes;
            treelet_boxes.reserve(treelets.size());
            for (auto const& el : treelets) {
                /*box3f box;
                box.lower = el.bounds_lower;
                box.upper = el.bounds_upper;*/
                treelet_boxes.push_back(el.bounds);
            }
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(
                bounds_data[pl_idx], treelet_boxes.data(), treelet_boxes.size() * sizeof(box3f), ctx.GetExecStream()));
        } else {
            CUDA_CHECK_ERROR(cuMemAllocAsync(&bounds_data[pl_idx], 1 * sizeof(box3f), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(
                cuMemcpyHtoDAsync(bounds_data[pl_idx], local_box.data(), 6 * sizeof(float), ctx.GetExecStream()));
        }


        //////////////////////////////////////
        // geometry
        //////////////////////////////////////


        buildInput.type = OPTIX_BUILD_INPUT_TYPE_CUSTOM_PRIMITIVES;
        auto& cp_input = buildInput.customPrimitiveArray;
        cp_input.aabbBuffers = &bounds_data[pl_idx];
        if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) {
            cp_input.numPrimitives = treelets.size();
        } else {
            cp_input.numPrimitives = 1;
        }
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
    ++geo_version;

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
    treelets_sbt_records_.clear();

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
        sbt_record.data.worldBounds.lower = std::get<0>(local_boxes_[pl_idx]);
        sbt_record.data.worldBounds.upper = std::get<1>(local_boxes_[pl_idx]);

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


        SBTRecord<device::TreeletsGeoData> treelets_sbt_record;
        OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(treelets_module_, &treelets_sbt_record));

        treelets_sbt_record.data.particleBufferPtr = (device::PKDParticle*) particle_data_[pl_idx];
        treelets_sbt_record.data.radiusBufferPtr = nullptr;
        treelets_sbt_record.data.colorBufferPtr = nullptr;
        treelets_sbt_record.data.treeletBufferPtr = (device::PKDlet*) treelets_data_[pl_idx];
        treelets_sbt_record.data.radius = particles.GetGlobalRadius();
        treelets_sbt_record.data.hasGlobalRadius = has_global_radius(particles);
        treelets_sbt_record.data.hasColorData = has_color(particles);
        treelets_sbt_record.data.globalColor =
            glm::vec4(particles.GetGlobalColour()[0] / 255.f, particles.GetGlobalColour()[1] / 255.f,
                particles.GetGlobalColour()[2] / 255.f, particles.GetGlobalColour()[3] / 255.f);
        treelets_sbt_record.data.particleCount = p_count;
        treelets_sbt_record.data.worldBounds.lower = std::get<0>(local_boxes_[pl_idx]);
        treelets_sbt_record.data.worldBounds.upper = std::get<1>(local_boxes_[pl_idx]);

        if (!has_global_radius(particles)) {
            treelets_sbt_record.data.radiusBufferPtr = (float*) radius_data_[pl_idx];
        }
        if (has_color(particles)) {
            treelets_sbt_record.data.colorBufferPtr = (glm::vec4*) color_data_[pl_idx];
        }
        treelets_sbt_records_.push_back(treelets_sbt_record);

        // occlusion stuff
        SBTRecord<device::TreeletsGeoData> treelets_sbt_record_occlusion;
        OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(treelets_occlusion_module_, &treelets_sbt_record_occlusion));

        treelets_sbt_record_occlusion.data = treelets_sbt_record.data;
        treelets_sbt_records_.push_back(treelets_sbt_record_occlusion);
    }

    ++sbt_version;

    return true;
}

} // namespace megamol::optix_hpg
