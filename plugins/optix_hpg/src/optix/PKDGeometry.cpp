#include "PKDGeometry.h"

#include <algorithm>
#include <fstream>
#include <mutex>
#include <unordered_set>

#include "mmcore/param/BoolParam.h"
#include "mmcore/param/EnumParam.h"
#ifndef MEGAMOL_USE_POWER
#include "mmcore/param/FilePathParam.h"
#endif
#include "mmcore/param/IntParam.h"

#include <glm/glm.hpp>

#include <omp.h>

#include <tbb/parallel_for.h>

#include "CallGeometry.h"

#include "PKDUtils.h"
#include "pkd_utils.h"

//#include "nvcomp.hpp"
//#include "nvcomp/lz4.hpp"
//#include "nvcomp/nvcompManagerFactory.hpp"

#include "SPKDGridify.h"

#include "QTreelets.h"

#include "moldyn/RDF.h"

#ifdef MEGAMOL_USE_POWER
#include <arrow/io/file.h>
#include <parquet/api/reader.h>
#include <parquet/api/writer.h>
#endif

namespace megamol::optix_hpg {
extern "C" const char embedded_pkd_programs[];

PKDGeometry::PKDGeometry()
        : out_geo_slot_("outGeo", "")
        , in_data_slot_("inData", "")
        , mode_slot_("mode", "")
        , compression_slot_("compression", "")
        , grid_slot_("grid", "")
        , threshold_slot_("threshold", "")
        , dump_debug_info_slot_("Debug::dumpDebugInfo", "")
#ifndef MEGAMOL_USE_POWER
        , debug_output_path_slot_("Debug::outputPath", "")
#endif
        , qtreelet_type_slot_("QTreelet::type", "") {
    out_geo_slot_.SetCallback(CallGeometry::ClassName(), CallGeometry::FunctionName(0), &PKDGeometry::get_data_cb);
    out_geo_slot_.SetCallback(CallGeometry::ClassName(), CallGeometry::FunctionName(1), &PKDGeometry::get_extents_cb);
    MakeSlotAvailable(&out_geo_slot_);

    in_data_slot_.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&in_data_slot_);

    auto ep = new core::param::EnumParam(static_cast<int>(PKDMode::STANDARD));
    ep->SetTypePair(static_cast<int>(PKDMode::STANDARD), "Standard");
    ep->SetTypePair(static_cast<int>(PKDMode::TREELETS), "Treelets");
    ep->SetTypePair(static_cast<int>(PKDMode::QTREELETS), "QTreelets");
    mode_slot_ << ep;
    MakeSlotAvailable(&mode_slot_);

    compression_slot_ << new core::param::BoolParam(false);
    MakeSlotAvailable(&compression_slot_);

    grid_slot_ << new core::param::BoolParam(false);
    MakeSlotAvailable(&grid_slot_);

    threshold_slot_ << new core::param::IntParam(256, 16);
    MakeSlotAvailable(&threshold_slot_);

    dump_debug_info_slot_ << new core::param::BoolParam(false);
    MakeSlotAvailable(&dump_debug_info_slot_);

#ifndef MEGAMOL_USE_POWER
    debug_output_path_slot_ << new core::param::FilePathParam("comp_stats.csv");
    MakeSlotAvailable(&debug_output_path_slot_);
#endif

    ep = new core::param::EnumParam(static_cast<int>(QTreeletType::E5M15D));
    ep->SetTypePair(static_cast<int>(QTreeletType::E5M15), "E5M15");
    ep->SetTypePair(static_cast<int>(QTreeletType::E4M16), "E4M16");
    ep->SetTypePair(static_cast<int>(QTreeletType::E4M16D), "E4M16D");
    ep->SetTypePair(static_cast<int>(QTreeletType::E5M15D), "E5M15D");
    qtreelet_type_slot_ << ep;
    MakeSlotAvailable(&qtreelet_type_slot_);
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
    /*for (auto& el : radius_data_) {
        CUDA_CHECK_ERROR(cuMemFree(el));
    }*/
    for (auto& el : color_data_) {
        CUDA_CHECK_ERROR(cuMemFree(el));
    }
    for (auto& el : treelets_data_) {
        CUDA_CHECK_ERROR(cuMemFree(el));
    }

    for (auto& el : exp_x_data_) {
        CUDA_CHECK_ERROR(cuMemFree(el));
    }
    for (auto& el : exp_y_data_) {
        CUDA_CHECK_ERROR(cuMemFree(el));
    }
    for (auto& el : exp_z_data_) {
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

    if (in_data->FrameID() != frame_id_ || in_data->DataHash() != data_hash_ || mode_slot_.IsDirty() ||
        ((mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) &&
            threshold_slot_.IsDirty()) ||
        ((mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) &&
            compression_slot_.IsDirty()) ||
        ((mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) &&
            compression_slot_.Param<core::param::BoolParam>()->Value() && grid_slot_.IsDirty()) ||
        ((mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::QTREELETS)) &&
            qtreelet_type_slot_.IsDirty())) {
        if (!assert_data(*in_data, *ctx))
            return false;
        createSBTRecords(*in_data, *ctx);
        frame_id_ = in_data->FrameID();
        data_hash_ = in_data->DataHash();
        if (mode_slot_.IsDirty() || threshold_slot_.IsDirty() ||
            ((mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) &&
                compression_slot_.IsDirty()) ||
            ((mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) &&
                compression_slot_.Param<core::param::BoolParam>()->Value() && grid_slot_.IsDirty()) ||
            ((mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::QTREELETS)) &&
                qtreelet_type_slot_.IsDirty())) {
            ++program_version;
        }
        mode_slot_.ResetDirty();
        threshold_slot_.ResetDirty();
        compression_slot_.ResetDirty();
        grid_slot_.ResetDirty();
        qtreelet_type_slot_.ResetDirty();
    }

    if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) {
        if (compression_slot_.Param<core::param::BoolParam>()->Value() &&
            !grid_slot_.Param<core::param::BoolParam>()->Value()) {
            program_groups_[0] = comp_treelets_module_;
            program_groups_[1] = comp_treelets_occlusion_module_;
        } else if (compression_slot_.Param<core::param::BoolParam>()->Value() &&
                   grid_slot_.Param<core::param::BoolParam>()->Value()) {
            // full compression
            program_groups_[0] = s_comp_treelets_module_;
            program_groups_[1] = s_comp_treelets_occlusion_module_;
        } else {
            program_groups_[0] = treelets_module_;
            program_groups_[1] = treelets_occlusion_module_;
        }
    } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::QTREELETS)) {
        program_groups_[0] = qpkd_treelets_module_;
        program_groups_[1] = qpkd_treelets_occlusion_module_;
    } else {
        program_groups_[0] = pkd_module_;
        program_groups_[1] = pkd_occlusion_module_;
    }

    out_geo->set_handle(&geo_handle_, geo_version);
    out_geo->set_program_groups(program_groups_.data(), program_groups_.size(), program_version);
    if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS)) {
        if (compression_slot_.Param<core::param::BoolParam>()->Value() &&
            !grid_slot_.Param<core::param::BoolParam>()->Value()) {
            out_geo->set_record(comp_treelets_sbt_records_.data(), comp_treelets_sbt_records_.size(),
                sizeof(SBTRecord<device::QTreeletsGeoData>), sbt_version);
        } else if (compression_slot_.Param<core::param::BoolParam>()->Value() &&
                   grid_slot_.Param<core::param::BoolParam>()->Value()) {
            // full compression
            out_geo->set_record(s_comp_treelets_sbt_records_.data(), s_comp_treelets_sbt_records_.size(),
                sizeof(SBTRecord<device::STreeletsGeoData>), sbt_version);
        } else {
            out_geo->set_record(treelets_sbt_records_.data(), treelets_sbt_records_.size(),
                sizeof(SBTRecord<device::TreeletsGeoData>), sbt_version);
        }
    } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::QTREELETS)) {
        out_geo->set_record(qpkd_treelets_sbt_records_.data(), qpkd_treelets_sbt_records_.size(),
            sizeof(SBTRecord<device::QPKDTreeletsGeoData>), sbt_version);
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

    pkd_occlusion_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(), &ctx.GetModuleCompileOptions(),
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

    comp_treelets_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(), &ctx.GetModuleCompileOptions(),
        &ctx.GetPipelineCompileOptions(), MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "comp_treelets_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "comp_treelets_closesthit"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});

    comp_treelets_occlusion_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
        &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
        MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "comp_treelets_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "comp_treelets_closesthit_occlusion"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});

    s_comp_treelets_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
        &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
        MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "s_comp_treelets_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "s_comp_treelets_closesthit"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});

    s_comp_treelets_occlusion_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
        &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
        MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
        {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "s_comp_treelets_intersect"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "s_comp_treelets_closesthit_occlusion"},
            {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});


    switch (static_cast<QTreeletType>(qtreelet_type_slot_.Param<core::param::EnumParam>()->Value())) {
    case QTreeletType::E5M15: {
        qpkd_treelets_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
            &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
            MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
            {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelet_intersect_e5m15"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "qpkd_treelets_closesthit"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});

        qpkd_treelets_occlusion_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
            &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
            MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
            {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelet_intersect_e5m15"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "qpkd_treelets_closesthit_occlusion"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});
    } break;
    case QTreeletType::E4M16: {
        qpkd_treelets_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
            &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
            MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
            {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelet_intersect_e4m16"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "qpkd_treelets_closesthit"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});

        qpkd_treelets_occlusion_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
            &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
            MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
            {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelet_intersect_e4m16"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "qpkd_treelets_closesthit_occlusion"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});
    } break;
    case QTreeletType::E5M15D: {
        qpkd_treelets_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
            &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
            MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
            {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelet_intersect_e5m15d"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "qpkd_treelets_closesthit"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});

        qpkd_treelets_occlusion_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
            &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
            MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
            {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelet_intersect_e5m15d"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "qpkd_treelets_closesthit_occlusion"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});
    } break;
    case QTreeletType::E4M16D: {
        qpkd_treelets_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
            &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
            MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
            {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelet_intersect_e4m16d"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "qpkd_treelets_closesthit"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});

        qpkd_treelets_occlusion_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
            &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
            MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
            {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelet_intersect_e4m16d"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "qpkd_treelets_closesthit_occlusion"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});
    } break;
    default: {
        qpkd_treelets_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
            &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
            MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
            {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelet_intersect_e4m16"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "qpkd_treelets_closesthit"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});

        qpkd_treelets_occlusion_module_ = MMOptixModule(embedded_pkd_programs, ctx.GetOptiXContext(),
            &ctx.GetModuleCompileOptions(), &ctx.GetPipelineCompileOptions(),
            MMOptixModule::MMOptixProgramGroupKind::MMOPTIX_PROGRAM_GROUP_KIND_HITGROUP,
            {{MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_INTERSECTION, "treelet_intersect_e4m16"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_CLOSESTHIT, "qpkd_treelets_closesthit_occlusion"},
                {MMOptixModule::MMOptixNameKind::MMOPTIX_NAME_BOUNDS, "treelets_bounds"}});
    } break;
    }


    ++program_version;

    return true;
}

void dump_analysis_data(std::filesystem::path const& output_path, std::shared_ptr<std::vector<glm::vec3>> op_s,
    std::shared_ptr<std::vector<glm::vec3>> sp_s, std::shared_ptr<std::vector<glm::vec3>> diffs, unsigned int pl_idx,
    float radius) {
    auto rdf = moldyn::RDF(op_s, sp_s);
    auto const [org_rdf, new_rdf] = rdf.BuildHistogram(4.0f * radius, 100);

    {
        auto f = std::ofstream(output_path / ("org_rdf_" + std::to_string(pl_idx) + ".blobb"));
        f.write(reinterpret_cast<char const*>(org_rdf.data()), org_rdf.size() * sizeof(decltype(org_rdf)::value_type));
        f.close();
    }

    {
        auto f = std::ofstream(output_path / ("new_rdf_" + std::to_string(pl_idx) + ".blobb"));
        f.write(reinterpret_cast<char const*>(new_rdf.data()), new_rdf.size() * sizeof(decltype(new_rdf)::value_type));
        f.close();
    }

    {
        auto const dx_minmax = std::minmax_element(
            diffs->begin(), diffs->end(), [](auto const& lhs, auto const& rhs) { return lhs.x < rhs.x; });
        auto const dy_minmax = std::minmax_element(
            diffs->begin(), diffs->end(), [](auto const& lhs, auto const& rhs) { return lhs.y < rhs.y; });
        auto const dz_minmax = std::minmax_element(
            diffs->begin(), diffs->end(), [](auto const& lhs, auto const& rhs) { return lhs.z < rhs.z; });
        auto const d_acc = std::accumulate(diffs->begin(), diffs->end(), glm::vec3(0));
        auto const csv_file_path = output_path / "comp_stats.csv";
        if (std::filesystem::exists(csv_file_path)) {
            // already exists ... append stats
            auto f = std::ofstream(csv_file_path, std::ios::app);
            f << dx_minmax.first->x << "," << dx_minmax.second->x << "," << dy_minmax.first->y << ","
              << dy_minmax.second->y << "," << dz_minmax.first->z << "," << dz_minmax.second->z << ","
              << d_acc.x / static_cast<float>(diffs->size()) << "," << d_acc.y / static_cast<float>(diffs->size())
              << "," << d_acc.z / static_cast<float>(diffs->size()) << "\n";
            f.close();
        } else {
            // create file
            auto f = std::ofstream(csv_file_path);
            f << "dx_min,dx_max,dy_min,dy_max,dz_min,dz_max,dx_mean,dy_mean,dz_mean\n";
            f << dx_minmax.first->x << "," << dx_minmax.second->x << "," << dy_minmax.first->y << ","
              << dy_minmax.second->y << "," << dz_minmax.first->z << "," << dz_minmax.second->z << ","
              << d_acc.x / static_cast<float>(diffs->size()) << "," << d_acc.y / static_cast<float>(diffs->size())
              << "," << d_acc.z / static_cast<float>(diffs->size()) << "\n";
            f.close();
        }
    }
}

bool PKDGeometry::assert_data(geocalls::MultiParticleDataCall const& call, Context const& ctx) {
#ifdef MEGAMOL_USE_POWER
    auto power_callbacks = this->frontend_resources.get<frontend_resources::PowerCallbacks>();
#endif

    free_data(ctx);

    auto const pl_count = call.GetParticleListCount();

    particle_data_.resize(pl_count, 0);
    //radius_data_.resize(pl_count, 0);
    color_data_.resize(pl_count, 0);
    treelets_data_.resize(pl_count, 0);
    local_boxes_.resize(pl_count);
    std::vector<CUdeviceptr> bounds_data(pl_count);
    std::vector<OptixBuildInput> build_inputs;

    exp_x_data_.resize(pl_count, 0);
    exp_y_data_.resize(pl_count, 0);
    exp_z_data_.resize(pl_count, 0);

    auto bbox = call.GetBoundingBoxes().ObjectSpaceBBox();
    glm::vec3 lower = glm::vec3(bbox.GetLeft(), bbox.Bottom(), bbox.Back());
    glm::vec3 upper = glm::vec3(bbox.GetRight(), bbox.Top(), bbox.Front());

    /*{
        std::ofstream coord_file = std::ofstream("./coord.csv");
        coord_file << "x,y,z,dx,dy,dz\n";
        coord_file.close();
    }*/

#ifdef MEGAMOL_USE_POWER
    size_t total_num_treelets = 0;
    size_t total_original_data_size = 0;
    size_t total_compressed_data_size = 0;
    size_t total_num_cells = 0;
#endif

    for (unsigned int pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
        auto const& particles = call.AccessParticles(pl_idx);

        auto const p_count = particles.GetCount();
        if (p_count == 0 || !has_global_radius(particles)) {
            if (!has_global_radius(particles)) {
                megamol::core::utility::log::Log::DefaultLog.WriteWarn(
                    "[PKDGeometry]: Per-particle radius not supported");
            }
            continue;
        }

        std::vector<device::PKDParticle> data(p_count);
        auto x_acc = particles.GetParticleStore().GetXAcc();
        auto y_acc = particles.GetParticleStore().GetYAcc();
        auto z_acc = particles.GetParticleStore().GetZAcc();

        // auto rad_acc = particles.GetParticleStore().GetRAcc();
        auto const global_rad = particles.GetGlobalRadius();

        auto cr_acc = particles.GetParticleStore().GetCRAcc();
        auto cg_acc = particles.GetParticleStore().GetCGAcc();
        auto cb_acc = particles.GetParticleStore().GetCBAcc();
        auto ca_acc = particles.GetParticleStore().GetCAAcc();

        for (std::size_t p_idx = 0; p_idx < p_count; ++p_idx) {
            data[p_idx].pos.x = x_acc->Get_f(p_idx);
            data[p_idx].pos.y = y_acc->Get_f(p_idx);
            data[p_idx].pos.z = z_acc->Get_f(p_idx);
        }

        device::box3f local_box;

        std::vector<device::PKDlet> treelets;
        std::vector<device::SPKDlet> s_treelets;
        std::vector<device::QPKDParticle> qparticles;
        std::vector<device::SPKDParticle> s_particles;

        if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS) &&
            compression_slot_.Param<core::param::BoolParam>()->Value() &&
            grid_slot_.Param<core::param::BoolParam>()->Value()) {
            // Grid Compression Treelets
            qparticles.resize(data.size());
            s_particles.resize(data.size());

            auto diffs = std::make_shared<std::vector<glm::vec3>>();
            diffs->reserve(data.size());
            auto orgpos = std::make_shared<std::vector<glm::vec3>>();
            orgpos->reserve(data.size());
            auto spos = std::make_shared<std::vector<glm::vec3>>();
            spos->reserve(data.size());

            /*= partition_data(
                data, threshold_slot_.Param<core::param::IntParam>()->Value(), particles.GetGlobalRadius());*/

            //convert(0, data.data(), qparticles.data(), data.size(), lower, global_rad);
            // separate into 256 grids
            //std::mutex data_add_mtx;
            auto const cells = gridify(data, lower, upper);
            megamol::core::utility::log::Log::DefaultLog.WriteInfo("[PKDGeometry] Num cells: %d", cells.size());
#ifdef MEGAMOL_USE_POWER
            total_num_cells += cells.size();
#endif
            for (auto const& c : cells) {
                auto const box = extendBounds(data, c.first, c.second, particles.GetGlobalRadius());
                auto tmp_t = partition_data(data, c.first, c.second, box.lower,
                    threshold_slot_.Param<core::param::IntParam>()->Value(), particles.GetGlobalRadius());
                for (auto& el : tmp_t) {
                    el.lower = box.lower;
                }
                /*std::transform(data.begin() + c.first, data.begin() + c.second, qparticles.begin() + c.first,
                    [&box](auto const& p) { return encode_coord(p.pos - box.lower, glm::vec3(), glm::vec3()); });*/
                //std::vector<device::SPKDParticle> tmp_p(c.second - c.first);
                for (auto const& el : tmp_t) {
                    std::transform(data.begin() + el.begin, data.begin() + el.end, s_particles.begin() + el.begin,
                        [&el, &box](auto const& p) {
                            auto const qp = encode_coord(p.pos - box.lower, glm::vec3(), glm::vec3());
                            device::SPKDParticle sp;
                            byte_cast bc;
                            bc.ui = 0;
                            bc.ui = qp.x;
                            sp.x = bc.parts.a;
                            auto fit_x = std::find(el.sx, el.sx + device::spkd_array_size, bc.parts.b);
                            if (fit_x == el.sx + device::spkd_array_size) {
                                throw std::runtime_error("did not find propper index");
                            }
                            sp.sx_idx = std::distance(el.sx, fit_x);
                            bc.ui = qp.y;
                            sp.y = bc.parts.a;
                            auto fit_y = std::find(el.sy, el.sy + device::spkd_array_size, bc.parts.b);
                            if (fit_y == el.sy + device::spkd_array_size) {
                                throw std::runtime_error("did not find propper index");
                            }
                            sp.sy_idx = std::distance(el.sy, fit_y);
                            bc.ui = qp.z;
                            sp.z = bc.parts.a;
                            auto fit_z = std::find(el.sz, el.sz + device::spkd_array_size, bc.parts.b);
                            if (fit_z == el.sz + device::spkd_array_size) {
                                throw std::runtime_error("did not find propper index");
                            }
                            sp.sz_idx = std::distance(el.sz, fit_z);
                            return sp;
                        });
                    //el.bounds = extendBounds(data, el.begin, el.end, particles.GetGlobalRadius());
                }

                if (dump_debug_info_slot_.Param<core::param::BoolParam>()->Value()) {
                    auto const [tmp_d, tmp_op, tmp_s] = compute_diffs(tmp_t, s_particles, data, c.first, c.second);
                    diffs->insert(diffs->end(), tmp_d.begin(), tmp_d.end());
                    orgpos->insert(orgpos->end(), tmp_op.begin(), tmp_op.end());
                    spos->insert(spos->end(), tmp_s.begin(), tmp_s.end());
                }

                // make PKD
                tbb::parallel_for(
                    (size_t) 0, tmp_t.size(), [&](size_t treeletID) { makePKD(s_particles, tmp_t[treeletID], 0); });
                s_treelets.insert(s_treelets.end(), tmp_t.begin(), tmp_t.end());
                //s_particles.insert(s_particles.end(), tmp_p.begin(), tmp_p.end());
            }

            CUDA_CHECK_ERROR(cuMemAllocAsync(
                &treelets_data_[pl_idx], s_treelets.size() * sizeof(device::SPKDlet), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(treelets_data_[pl_idx], s_treelets.data(),
                s_treelets.size() * sizeof(device::SPKDlet), ctx.GetExecStream()));
            megamol::core::utility::log::Log::DefaultLog.WriteInfo("[PKDGeometry] Num treelets: %d", s_treelets.size());
            megamol::core::utility::log::Log::DefaultLog.WriteInfo("[PKDGeometry] Original size: %d; New size: %d",
                data.size() * sizeof(device::PKDParticle),
                s_treelets.size() * sizeof(device::SPKDlet) + s_particles.size() * sizeof(device::SPKDParticle));

#ifdef MEGAMOL_USE_POWER
            total_num_treelets += s_treelets.size();
            total_original_data_size += data.size() * sizeof(glm::vec3);
            total_compressed_data_size +=
                s_treelets.size() * sizeof(device::SPKDlet) + s_particles.size() * sizeof(device::SPKDParticle);
#endif

            if (dump_debug_info_slot_.Param<core::param::BoolParam>()->Value()) {
#ifdef MEGAMOL_USE_POWER
                auto const output_path = power_callbacks.get_output_path();
#else
                auto const output_path = debug_output_path_slot_.Param<core::param::FilePathParam>()->Value();
#endif
                dump_analysis_data(output_path, orgpos, spos, diffs, pl_idx, particles.GetGlobalRadius());
            }

#if 0
                auto op_s = std::make_shared<std::vector<glm::vec3>>(orgpos.begin(), orgpos.end());
                auto sp_s = std::make_shared<std::vector<glm::vec3>>(spos.begin(), spos.end());
                auto rdf = moldyn::RDF(op_s, sp_s);
                auto const [org_rdf, new_rdf] = rdf.BuildHistogram(4.0f * particles.GetGlobalRadius(), 100);

                {
                    auto f = std::ofstream(output_path / ("org_rdf_" + std::to_string(pl_idx) + ".blobb"));
                    f.write(reinterpret_cast<char const*>(org_rdf.data()),
                        org_rdf.size() * sizeof(decltype(org_rdf)::value_type));
                    f.close();
                }

                {
                    auto f = std::ofstream(output_path / ("new_rdf_" + std::to_string(pl_idx) + ".blobb"));
                    f.write(reinterpret_cast<char const*>(new_rdf.data()),
                        new_rdf.size() * sizeof(decltype(new_rdf)::value_type));
                    f.close();
                }

                {
                    auto const dx_minmax = std::minmax_element(
                        diffs.begin(), diffs.end(), [](auto const& lhs, auto const& rhs) { return lhs.x < rhs.x; });
                    auto const dy_minmax = std::minmax_element(
                        diffs.begin(), diffs.end(), [](auto const& lhs, auto const& rhs) { return lhs.y < rhs.y; });
                    auto const dz_minmax = std::minmax_element(
                        diffs.begin(), diffs.end(), [](auto const& lhs, auto const& rhs) { return lhs.z < rhs.z; });
                    auto const d_acc = std::accumulate(diffs.begin(), diffs.end(), glm::vec3(0));
                    auto const csv_file_path = output_path / "comp_stats.csv";
                    if (std::filesystem::exists(csv_file_path)) {
                        // already exists ... append stats
                        auto f = std::ofstream(csv_file_path, std::ios::app);
                        f << dx_minmax.first->x << "," << dx_minmax.second->x << "," << dy_minmax.first->y << ","
                          << dy_minmax.second->y << "," << dz_minmax.first->z << "," << dz_minmax.second->z << ","
                          << d_acc.x / static_cast<float>(diffs.size()) << ","
                          << d_acc.y / static_cast<float>(diffs.size()) << ","
                          << d_acc.z / static_cast<float>(diffs.size()) << "\n";
                        f.close();
                    } else {
                        // create file
                        auto f = std::ofstream(csv_file_path);
                        f << "dx_min,dx_max,dy_min,dy_max,dz_min,dz_max,dx_mean,dy_mean,dz_mean\n";
                        f << dx_minmax.first->x << "," << dx_minmax.second->x << "," << dy_minmax.first->y << ","
                          << dy_minmax.second->y << "," << dz_minmax.first->z << "," << dz_minmax.second->z << ","
                          << d_acc.x / static_cast<float>(diffs.size()) << ","
                          << d_acc.y / static_cast<float>(diffs.size()) << ","
                          << d_acc.z / static_cast<float>(diffs.size()) << "\n";
                        f.close();
                    }
                }
#endif

#if 0
#ifdef MEGAMOL_USE_POWER
                auto const file_path = output_path / "diff.parquet";
                {
                    using namespace parquet;
                    using namespace parquet::schema;

                    try {
                        std::size_t min_field_size = 0;
                        bool first_time = true;

                        // create scheme
                        NodeVector fields;
                        fields.reserve(9);
                        fields.push_back(
                            PrimitiveNode::Make("x", Repetition::REQUIRED, Type::FLOAT, ConvertedType::NONE));
                        fields.push_back(
                            PrimitiveNode::Make("y", Repetition::REQUIRED, Type::FLOAT, ConvertedType::NONE));
                        fields.push_back(
                            PrimitiveNode::Make("z", Repetition::REQUIRED, Type::FLOAT, ConvertedType::NONE));
                        fields.push_back(
                            PrimitiveNode::Make("sx", Repetition::REQUIRED, Type::FLOAT, ConvertedType::NONE));
                        fields.push_back(
                            PrimitiveNode::Make("sy", Repetition::REQUIRED, Type::FLOAT, ConvertedType::NONE));
                        fields.push_back(
                            PrimitiveNode::Make("sz", Repetition::REQUIRED, Type::FLOAT, ConvertedType::NONE));
                        fields.push_back(
                            PrimitiveNode::Make("dx", Repetition::REQUIRED, Type::FLOAT, ConvertedType::NONE));
                        fields.push_back(
                            PrimitiveNode::Make("dy", Repetition::REQUIRED, Type::FLOAT, ConvertedType::NONE));
                        fields.push_back(
                            PrimitiveNode::Make("dz", Repetition::REQUIRED, Type::FLOAT, ConvertedType::NONE));
                        auto schema = std::static_pointer_cast<GroupNode>(
                            GroupNode::Make("schema", Repetition::REQUIRED, fields));

                        // open file
                        std::shared_ptr<::arrow::io::FileOutputStream> file;
                        PARQUET_ASSIGN_OR_THROW(file, ::arrow::io::FileOutputStream::Open(file_path.string()));

                        // configure
                        WriterProperties::Builder builder;
                        builder.compression(Compression::BROTLI);
                        auto props = builder.build();

                        // create instance
                        auto file_writer = ParquetFileWriter::Open(file, schema, props);

                        // write data
                        auto rg_writer = file_writer->AppendBufferedRowGroup();

                        std::vector<float> x_vals(orgpos.size());
                        std::transform(orgpos.begin(), orgpos.end(), x_vals.begin(), [](auto const& d) { return d.x; });
                        std::vector<float> y_vals(orgpos.size());
                        std::transform(orgpos.begin(), orgpos.end(), y_vals.begin(), [](auto const& d) { return d.y; });
                        std::vector<float> z_vals(orgpos.size());
                        std::transform(orgpos.begin(), orgpos.end(), z_vals.begin(), [](auto const& d) { return d.z; });

                        auto float_writer = static_cast<parquet::FloatWriter*>(rg_writer->column(0));
                        float_writer->WriteBatch(x_vals.size(), nullptr, nullptr, x_vals.data());
                        float_writer = static_cast<parquet::FloatWriter*>(rg_writer->column(1));
                        float_writer->WriteBatch(y_vals.size(), nullptr, nullptr, y_vals.data());
                        float_writer = static_cast<parquet::FloatWriter*>(rg_writer->column(2));
                        float_writer->WriteBatch(z_vals.size(), nullptr, nullptr, z_vals.data());

                        x_vals.resize(spos.size());
                        std::transform(spos.begin(), spos.end(), x_vals.begin(), [](auto const& d) { return d.x; });
                        y_vals.resize(spos.size());
                        std::transform(spos.begin(), spos.end(), y_vals.begin(), [](auto const& d) { return d.y; });
                        z_vals.resize(spos.size());
                        std::transform(spos.begin(), spos.end(), z_vals.begin(), [](auto const& d) { return d.z; });

                        float_writer = static_cast<parquet::FloatWriter*>(rg_writer->column(3));
                        float_writer->WriteBatch(x_vals.size(), nullptr, nullptr, x_vals.data());
                        float_writer = static_cast<parquet::FloatWriter*>(rg_writer->column(4));
                        float_writer->WriteBatch(y_vals.size(), nullptr, nullptr, y_vals.data());
                        float_writer = static_cast<parquet::FloatWriter*>(rg_writer->column(5));
                        float_writer->WriteBatch(z_vals.size(), nullptr, nullptr, z_vals.data());

                        x_vals.resize(diffs.size());
                        std::transform(diffs.begin(), diffs.end(), x_vals.begin(), [](auto const& d) { return d.x; });
                        y_vals.resize(diffs.size());
                        std::transform(diffs.begin(), diffs.end(), y_vals.begin(), [](auto const& d) { return d.y; });
                        z_vals.resize(diffs.size());
                        std::transform(diffs.begin(), diffs.end(), z_vals.begin(), [](auto const& d) { return d.z; });

                        float_writer = static_cast<parquet::FloatWriter*>(rg_writer->column(6));
                        float_writer->WriteBatch(x_vals.size(), nullptr, nullptr, x_vals.data());
                        float_writer = static_cast<parquet::FloatWriter*>(rg_writer->column(7));
                        float_writer->WriteBatch(y_vals.size(), nullptr, nullptr, y_vals.data());
                        float_writer = static_cast<parquet::FloatWriter*>(rg_writer->column(8));
                        float_writer->WriteBatch(z_vals.size(), nullptr, nullptr, z_vals.data());

                        // close
                        rg_writer->Close();
                        file_writer->Close();
                    } catch (std::exception const& ex) {
                        core::utility::log::Log::DefaultLog.WriteError("[ParquetWriter]: %s", ex.what());
                    }
                }
#endif
            }
#endif
        }

        std::shared_ptr<QTPBufferBase> qtpbuffer;
        std::vector<device::QPKDlet> qtreelets;
        std::vector<char> exp_vec_x;
        std::vector<char> exp_vec_y;
        std::vector<char> exp_vec_z;
        if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::QTREELETS)) {
            // QTreelets
            auto const selected_type =
                static_cast<QTreeletType>(qtreelet_type_slot_.Param<core::param::EnumParam>()->Value());

            // 1 partion particles
            treelets = prePartition_inPlace(
                data, threshold_slot_.Param<core::param::IntParam>()->Value(), particles.GetGlobalRadius());

            // 2 make PKDs and quantize
            qtreelets.resize(treelets.size());
            std::vector<device::QTParticle> qtparticles(data.size());

            std::vector<unsigned int> bucket_count(treelets.size());
            std::vector<unsigned int> global_histo(256, 0);
            std::mutex histo_guard;

            tbb::parallel_for((size_t) 0, treelets.size(), [&](size_t treeletID) {
                makePKD(data, treelets[treeletID].begin, treelets[treeletID].end, treelets[treeletID].bounds);

                auto const histo_t = quantizeTree2(
                    selected_type, data.data(), treelets[treeletID], qtparticles.data(), qtreelets[treeletID]);

                {
                    auto num_el = std::count_if(std::get<0>(histo_t).begin(), std::get<0>(histo_t).end(),
                        [](auto const& val) { return val != 0; });
                    num_el = std::max(num_el, std::count_if(std::get<1>(histo_t).begin(), std::get<1>(histo_t).end(),
                                                  [](auto const& val) { return val != 0; }));
                    num_el = std::max(num_el, std::count_if(std::get<2>(histo_t).begin(), std::get<2>(histo_t).end(),
                                                  [](auto const& val) { return val != 0; }));
                    bucket_count[treeletID] = num_el;
                    std::lock_guard<std::mutex> histo_lock(histo_guard);
                    std::transform(std::get<0>(histo_t).begin(), std::get<0>(histo_t).end(), global_histo.begin(),
                        global_histo.begin(), [](auto vala, auto valb) { return vala + valb; });
                    std::transform(std::get<1>(histo_t).begin(), std::get<1>(histo_t).end(), global_histo.begin(),
                        global_histo.begin(), [](auto vala, auto valb) { return vala + valb; });
                    std::transform(std::get<2>(histo_t).begin(), std::get<2>(histo_t).end(), global_histo.begin(),
                        global_histo.begin(), [](auto vala, auto valb) { return vala + valb; });
                }
            });

            // 3 create exponent maps
            std::unordered_set<char> exponents_x;
            std::unordered_set<char> exponents_y;
            std::unordered_set<char> exponents_z;

            unsigned num_idx = 0;
            switch (selected_type) {
            case QTreeletType::E5M15:
                num_idx = static_cast<unsigned int>(std::pow(2, device::QTParticle_e5m15::exp));
                break;
            case QTreeletType::E4M16:
                num_idx = static_cast<unsigned int>(std::pow(2, device::QTParticle_e4m16::exp));
                break;
            case QTreeletType::E5M15D:
                num_idx = static_cast<unsigned int>(std::pow(2, device::QTParticle_e5m15d::exp));
                break;
            case QTreeletType::E4M16D:
                num_idx = static_cast<unsigned int>(std::pow(2, device::QTParticle_e4m16d::exp));
                break;
            default:
                num_idx = static_cast<unsigned int>(std::pow(2, device::QTParticle_e4m16::exp));
                break;
            }

            create_exp_maps(qtparticles, exponents_x, exponents_y, exponents_z, num_idx);

            exp_vec_x.resize(treelets.size() * num_idx);
            exp_vec_y.resize(treelets.size() * num_idx);
            exp_vec_z.resize(treelets.size() * num_idx);
            exp_vec_x = std::vector<char>(exponents_x.begin(), exponents_x.end());
            exp_vec_y = std::vector<char>(exponents_y.begin(), exponents_y.end());
            exp_vec_z = std::vector<char>(exponents_z.begin(), exponents_z.end());

            // 4 convert to qlets
            switch (selected_type) {
            case QTreeletType::E5M15:
                qtpbuffer = std::make_shared<QTPBuffer_e5m15>(data.size());
                break;
            case QTreeletType::E4M16:
                qtpbuffer = std::make_shared<QTPBuffer_e4m16>(data.size());
                break;
            case QTreeletType::E5M15D:
                qtpbuffer = std::make_shared<QTPBuffer_e5m15d>(data.size());
                break;
            case QTreeletType::E4M16D:
                qtpbuffer = std::make_shared<QTPBuffer_e4m16d>(data.size());
                break;
            default:
                qtpbuffer = std::make_shared<QTPBuffer_e4m16>(data.size());
            }

            switch (selected_type) {
            case QTreeletType::E5M15: {
                tbb::parallel_for((size_t) 0, treelets.size(), [&](size_t treeletID) {
                    unsigned int offset = 0;
                    convert_qlet<device::QTParticle_e5m15>(qtreelets[treeletID], qtparticles,
                        std::dynamic_pointer_cast<QTPBuffer_e5m15>(qtpbuffer)->buffer, exp_vec_x.data() + offset,
                        exp_vec_y.data() + offset, exp_vec_z.data() + offset);
                });
            } break;
            case QTreeletType::E4M16: {
                tbb::parallel_for((size_t) 0, treelets.size(), [&](size_t treeletID) {
                    unsigned int offset = 0;
                    convert_qlet<device::QTParticle_e4m16>(qtreelets[treeletID], qtparticles,
                        std::dynamic_pointer_cast<QTPBuffer_e4m16>(qtpbuffer)->buffer, exp_vec_x.data() + offset,
                        exp_vec_y.data() + offset, exp_vec_z.data() + offset);
                });
            } break;
            case QTreeletType::E5M15D: {
                tbb::parallel_for((size_t) 0, treelets.size(), [&](size_t treeletID) {
                    unsigned int offset = 0;
                    convert_qlet_dep<device::QTParticle_e5m15d>(qtreelets[treeletID], qtparticles,
                        std::dynamic_pointer_cast<QTPBuffer_e5m15d>(qtpbuffer)->buffer, exp_vec_x.data() + offset,
                        exp_vec_y.data() + offset, exp_vec_z.data() + offset);
                });
            } break;
            case QTreeletType::E4M16D: {
                tbb::parallel_for((size_t) 0, treelets.size(), [&](size_t treeletID) {
                    unsigned int offset = 0;
                    convert_qlet_dep<device::QTParticle_e4m16d>(qtreelets[treeletID], qtparticles,
                        std::dynamic_pointer_cast<QTPBuffer_e4m16d>(qtpbuffer)->buffer, exp_vec_x.data() + offset,
                        exp_vec_y.data() + offset, exp_vec_z.data() + offset);
                });
            } break;
            default:
                std::cout << "Should not happen" << std::endl;
            }

            CUDA_CHECK_ERROR(cuMemAllocAsync(
                &treelets_data_[pl_idx], qtreelets.size() * sizeof(device::QPKDlet), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(treelets_data_[pl_idx], qtreelets.data(),
                qtreelets.size() * sizeof(device::QPKDlet), ctx.GetExecStream()));

            if (dump_debug_info_slot_.Param<core::param::BoolParam>()->Value()) {
#ifdef MEGAMOL_USE_POWER
                auto const output_path = power_callbacks.get_output_path();
#else
                auto const output_path = debug_output_path_slot_.Param<core::param::FilePathParam>()->Value();
#endif
                auto diffs = std::make_shared<std::vector<glm::vec3>>();
                diffs->reserve(data.size());
                auto orgpos = std::make_shared<std::vector<glm::vec3>>();
                orgpos->reserve(data.size());
                auto newpos = std::make_shared<std::vector<glm::vec3>>();
                newpos->reserve(data.size());

                for (size_t i = 0; i < treelets.size(); ++i) {
                    unsigned int offset = 0;

                    /*if (local_tables) {
                        offset = i * num_idx;
                    }*/

                    auto const [diffs_t, orgpos_t, newpos_t] =
                        unified_sub_print(selected_type, 0, qtreelets[i].basePos, data.data(), qtpbuffer, qtreelets[i],
                            exp_vec_x.data() + offset, exp_vec_y.data() + offset, exp_vec_z.data() + offset);

                    diffs->insert(diffs->end(), diffs_t.begin(), diffs_t.end());
                    orgpos->insert(orgpos->end(), orgpos_t.begin(), orgpos_t.end());
                    newpos->insert(newpos->end(), newpos_t.begin(), newpos_t.end());
                }

                dump_analysis_data(output_path, orgpos, newpos, diffs, pl_idx, particles.GetGlobalRadius());
            }

#ifdef MEGAMOL_USE_POWER
            total_num_treelets += qtreelets.size();
            total_original_data_size += data.size() * sizeof(glm::vec3);
            total_compressed_data_size +=
                qtreelets.size() * sizeof(device::QPKDlet) + data.size() * sizeof(device::QTParticle_e5m15);
#endif
        }

        if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS) &&
            !(compression_slot_.Param<core::param::BoolParam>()->Value() &&
                grid_slot_.Param<core::param::BoolParam>()->Value())) {
            // Normal treelets
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
            megamol::core::utility::log::Log::DefaultLog.WriteInfo("[PKDGeometry] Num treelets: %d", treelets.size());

#ifdef MEGAMOL_USE_POWER
            total_num_treelets += treelets.size();
            total_original_data_size += data.size() * sizeof(glm::vec3);
            total_compressed_data_size +=
                treelets.size() * sizeof(device::PKDlet) + data.size() * sizeof(device::PKDParticle);
            /*power_callbacks.add_meta_key_value("NumTreelets", std::to_string(treelets.size()));
            power_callbacks.add_meta_key_value(
                "OriginalDataSize", std::to_string(data.size() * sizeof(device::PKDParticle)));
            power_callbacks.add_meta_key_value("CompressedDataSize",
                std::to_string(treelets.size() * sizeof(device::PKDlet) + data.size() * sizeof(device::PKDParticle)));*/
#endif

            // TODO compress data if requested
            // for debugging without parallel
            if (compression_slot_.Param<core::param::BoolParam>()->Value() &&
                !grid_slot_.Param<core::param::BoolParam>()->Value()) {
                // Normal treelets with compression
                qparticles.resize(data.size());
                for (size_t tID = 0; tID < treelets.size(); ++tID) {
                    auto const& treelet = treelets[tID];
                    std::vector<glm::uvec3> out_coord(treelet.end - treelet.begin);
                    std::vector<device::PKDParticle> out_decode(treelet.end - treelet.begin);
                    convert(0, &data[treelet.begin], &qparticles[treelet.begin], treelet.end - treelet.begin,
                        treelet.bounds, global_rad, out_decode.data(), out_coord.data());
                }
            }
        } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::STANDARD)) {
            // Normal PKDs
            auto const max_threads = omp_get_max_threads();
            std::vector<device::box3f> local_boxes(max_threads);
#pragma omp parallel for shared(local_boxes)
            for (int64_t p_idx = 0; p_idx < p_count; ++p_idx) {
                auto const thread_num = omp_get_thread_num();
                auto& box = local_boxes[thread_num];
                glm::vec3 pos(x_acc->Get_f(p_idx), y_acc->Get_f(p_idx), z_acc->Get_f(p_idx));
                auto const new_lower = pos - global_rad;
                auto const new_upper = pos + global_rad;
                box.extend(new_lower);
                box.extend(new_upper);
            }

            for (auto const& el : local_boxes) {
                local_box.extend(el);
            }

            makePKD(data, local_box);
            local_boxes_[pl_idx] = local_box;
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
        if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS) &&
            (compression_slot_.Param<core::param::BoolParam>()->Value() &&
                !grid_slot_.Param<core::param::BoolParam>()->Value())) {
            CUDA_CHECK_ERROR(
                cuMemAllocAsync(&particle_data_[pl_idx], p_count * sizeof(device::QPKDParticle), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(particle_data_[pl_idx], qparticles.data(),
                p_count * sizeof(device::QPKDParticle), ctx.GetExecStream()));
        } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS) &&
                   (compression_slot_.Param<core::param::BoolParam>()->Value() &&
                       grid_slot_.Param<core::param::BoolParam>()->Value())) {
            CUDA_CHECK_ERROR(
                cuMemAllocAsync(&particle_data_[pl_idx], p_count * sizeof(device::SPKDParticle), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(particle_data_[pl_idx], s_particles.data(),
                p_count * sizeof(device::SPKDParticle), ctx.GetExecStream()));
        } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::QTREELETS)) {
            CUDA_CHECK_ERROR(cuMemAllocAsync(&exp_x_data_[pl_idx], exp_vec_x.size(), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(cuMemAllocAsync(&exp_y_data_[pl_idx], exp_vec_y.size(), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(cuMemAllocAsync(&exp_z_data_[pl_idx], exp_vec_z.size(), ctx.GetExecStream()));

            CUDA_CHECK_ERROR(
                cuMemcpyHtoDAsync(exp_x_data_[pl_idx], exp_vec_x.data(), exp_vec_x.size(), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(
                cuMemcpyHtoDAsync(exp_y_data_[pl_idx], exp_vec_y.data(), exp_vec_y.size(), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(
                cuMemcpyHtoDAsync(exp_z_data_[pl_idx], exp_vec_z.data(), exp_vec_z.size(), ctx.GetExecStream()));

            switch (static_cast<QTreeletType>(qtreelet_type_slot_.Param<core::param::EnumParam>()->Value())) {
            case QTreeletType::E5M15: {
                auto const buf = std::dynamic_pointer_cast<QTPBuffer_e5m15>(qtpbuffer);
                CUDA_CHECK_ERROR(cuMemAllocAsync(
                    &particle_data_[pl_idx], p_count * sizeof(device::QTParticle_e5m15), ctx.GetExecStream()));
                CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(particle_data_[pl_idx], buf->buffer.data(),
                    buf->buffer.size() * sizeof(device::QTParticle_e5m15), ctx.GetExecStream()));
                /*total_size += buf->buffer.size() * sizeof(buf->buffer[0]);
                total_size_wo_treelets += buf->buffer.size() * sizeof(buf->buffer[0]);*/
            } break;
            case QTreeletType::E4M16: {
                auto const buf = std::dynamic_pointer_cast<QTPBuffer_e4m16>(qtpbuffer);
                CUDA_CHECK_ERROR(cuMemAllocAsync(
                    &particle_data_[pl_idx], p_count * sizeof(device::QTParticle_e4m16), ctx.GetExecStream()));
                CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(particle_data_[pl_idx], buf->buffer.data(),
                    buf->buffer.size() * sizeof(device::QTParticle_e4m16), ctx.GetExecStream()));
                /*total_size += buf->buffer.size() * sizeof(buf->buffer[0]);
                total_size_wo_treelets += buf->buffer.size() * sizeof(buf->buffer[0]);*/
            } break;
            case QTreeletType::E5M15D: {
                auto const buf = std::dynamic_pointer_cast<QTPBuffer_e5m15d>(qtpbuffer);
                CUDA_CHECK_ERROR(cuMemAllocAsync(
                    &particle_data_[pl_idx], p_count * sizeof(device::QTParticle_e5m15d), ctx.GetExecStream()));
                CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(particle_data_[pl_idx], buf->buffer.data(),
                    buf->buffer.size() * sizeof(device::QTParticle_e5m15d), ctx.GetExecStream()));
                /*total_size += buf->buffer.size() * sizeof(buf->buffer[0]);
                total_size_wo_treelets += buf->buffer.size() * sizeof(buf->buffer[0]);*/
            } break;
            case QTreeletType::E4M16D: {
                auto const buf = std::dynamic_pointer_cast<QTPBuffer_e4m16d>(qtpbuffer);
                CUDA_CHECK_ERROR(cuMemAllocAsync(
                    &particle_data_[pl_idx], p_count * sizeof(device::QTParticle_e4m16d), ctx.GetExecStream()));
                CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(particle_data_[pl_idx], buf->buffer.data(),
                    buf->buffer.size() * sizeof(device::QTParticle_e4m16d), ctx.GetExecStream()));
                /*total_size += buf->buffer.size() * sizeof(buf->buffer[0]);
                total_size_wo_treelets += buf->buffer.size() * sizeof(buf->buffer[0]);*/
            } break;
            default:
                std::cout << "Should not happen" << std::endl;
            }
        } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::STANDARD)) {
            CUDA_CHECK_ERROR(
                cuMemAllocAsync(&particle_data_[pl_idx], p_count * sizeof(device::PKDParticle), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(
                particle_data_[pl_idx], data.data(), p_count * sizeof(device::PKDParticle), ctx.GetExecStream()));
        }

        /*CUDA_CHECK_ERROR(cuMemAllocAsync(&radius_data_[pl_idx], rad_data.size() * sizeof(float), ctx.GetExecStream()));
        CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(
            radius_data_[pl_idx], rad_data.data(), rad_data.size() * sizeof(float), ctx.GetExecStream()));*/

        unsigned int geo_flag = OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT;
        build_inputs.emplace_back();
        OptixBuildInput& buildInput = build_inputs.back();
        memset(&buildInput, 0, sizeof(OptixBuildInput));

        if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS) &&
            !grid_slot_.Param<core::param::BoolParam>()->Value()) {
            // TODO set of treelet boxes
            CUDA_CHECK_ERROR(
                cuMemAllocAsync(&bounds_data[pl_idx], treelets.size() * sizeof(device::box3f), ctx.GetExecStream()));
            std::vector<device::box3f> treelet_boxes;
            treelet_boxes.reserve(treelets.size());
            for (auto const& el : treelets) {
                /*box3f box;
                box.lower = el.bounds_lower;
                box.upper = el.bounds_upper;*/
                treelet_boxes.push_back(el.bounds);
            }
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(bounds_data[pl_idx], treelet_boxes.data(),
                treelet_boxes.size() * sizeof(device::box3f), ctx.GetExecStream()));
        } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS) &&
                   grid_slot_.Param<core::param::BoolParam>()->Value()) {
            CUDA_CHECK_ERROR(
                cuMemAllocAsync(&bounds_data[pl_idx], s_treelets.size() * sizeof(device::box3f), ctx.GetExecStream()));
            std::vector<device::box3f> treelet_boxes;
            treelet_boxes.reserve(treelets.size());
            for (auto const& el : s_treelets) {
                /*box3f box;
                box.lower = el.bounds_lower;
                box.upper = el.bounds_upper;*/
                treelet_boxes.push_back(el.bounds);
            }
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(bounds_data[pl_idx], treelet_boxes.data(),
                treelet_boxes.size() * sizeof(device::box3f), ctx.GetExecStream()));
        } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::QTREELETS)) {
            CUDA_CHECK_ERROR(
                cuMemAllocAsync(&bounds_data[pl_idx], qtreelets.size() * sizeof(device::box3f), ctx.GetExecStream()));
            std::vector<device::box3f> treelet_boxes;
            treelet_boxes.reserve(qtreelets.size());
            for (auto const& el : qtreelets) {
                /*box3f box;
                box.lower = el.bounds_lower;
                box.upper = el.bounds_upper;*/
                treelet_boxes.push_back(el.bounds);
            }
            CUDA_CHECK_ERROR(cuMemcpyHtoDAsync(bounds_data[pl_idx], treelet_boxes.data(),
                treelet_boxes.size() * sizeof(device::box3f), ctx.GetExecStream()));
        } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::STANDARD)) {
            CUDA_CHECK_ERROR(cuMemAllocAsync(&bounds_data[pl_idx], 1 * sizeof(device::box3f), ctx.GetExecStream()));
            CUDA_CHECK_ERROR(
                cuMemcpyHtoDAsync(bounds_data[pl_idx], &local_box, sizeof(local_box), ctx.GetExecStream()));
        }


        //////////////////////////////////////
        // geometry
        //////////////////////////////////////


        buildInput.type = OPTIX_BUILD_INPUT_TYPE_CUSTOM_PRIMITIVES;
        auto& cp_input = buildInput.customPrimitiveArray;
        cp_input.aabbBuffers = &bounds_data[pl_idx];
        if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS) &&
            !grid_slot_.Param<core::param::BoolParam>()->Value()) {
            cp_input.numPrimitives = treelets.size();
        } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::TREELETS) &&
                   grid_slot_.Param<core::param::BoolParam>()->Value()) {
            cp_input.numPrimitives = s_treelets.size();
        } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(PKDMode::QTREELETS)) {
            cp_input.numPrimitives = qtreelets.size();
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
    CUDA_CHECK_ERROR(cuMemFreeAsync(geo_temp, ctx.GetExecStream()));

    std::size_t compSize;
    CUDA_CHECK_ERROR(cuMemcpyDtoHAsync(&compSize, d_compSize, sizeof(std::size_t), ctx.GetExecStream()));
    CUDA_CHECK_ERROR(cuMemFreeAsync(d_compSize, ctx.GetExecStream()));
#ifdef MEGAMOL_USE_POWER
    if (compSize < bufferSizes.outputSizeInBytes) {
        power_callbacks.add_meta_key_value("GeoSize", std::to_string(compSize));
    } else {
        power_callbacks.add_meta_key_value("GeoSize", std::to_string(bufferSizes.outputSizeInBytes));
    }
    power_callbacks.add_meta_key_value("NumTreelets", std::to_string(total_num_treelets));
    power_callbacks.add_meta_key_value("OriginalDataSize", std::to_string(total_original_data_size));
    power_callbacks.add_meta_key_value("CompressedDataSize", std::to_string(total_compressed_data_size));
    power_callbacks.add_meta_key_value("NumGridCells", std::to_string(total_num_cells));
#endif
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
    comp_treelets_sbt_records_.clear();
    s_comp_treelets_sbt_records_.clear();
    qpkd_treelets_sbt_records_.clear();

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
        sbt_record.data.worldBounds = local_boxes_[pl_idx];

        /*if (!has_global_radius(particles)) {
            sbt_record.data.radiusBufferPtr = (float*) radius_data_[pl_idx];
        }*/
        sbt_record.data.radiusBufferPtr = nullptr;
        if (has_color(particles)) {
            sbt_record.data.colorBufferPtr = (glm::vec4*) color_data_[pl_idx];
        }
        sbt_records_.push_back(sbt_record);

        // occlusion stuff
        SBTRecord<device::PKDGeoData> sbt_record_occlusion;
        OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(pkd_occlusion_module_, &sbt_record_occlusion));

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
        treelets_sbt_record.data.worldBounds = local_boxes_[pl_idx];

        /*if (!has_global_radius(particles)) {
            treelets_sbt_record.data.radiusBufferPtr = (float*) radius_data_[pl_idx];
        }*/
        treelets_sbt_record.data.radiusBufferPtr = nullptr;
        if (has_color(particles)) {
            treelets_sbt_record.data.colorBufferPtr = (glm::vec4*) color_data_[pl_idx];
        }
        treelets_sbt_records_.push_back(treelets_sbt_record);

        // occlusion stuff
        SBTRecord<device::TreeletsGeoData> treelets_sbt_record_occlusion;
        OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(treelets_occlusion_module_, &treelets_sbt_record_occlusion));

        treelets_sbt_record_occlusion.data = treelets_sbt_record.data;
        treelets_sbt_records_.push_back(treelets_sbt_record_occlusion);


        SBTRecord<device::QTreeletsGeoData> comp_treelets_sbt_record;
        OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(comp_treelets_module_, &comp_treelets_sbt_record));

        comp_treelets_sbt_record.data.particleBufferPtr = (device::QPKDParticle*) particle_data_[pl_idx];
        //comp_treelets_sbt_record.data.radiusBufferPtr = nullptr;
        comp_treelets_sbt_record.data.colorBufferPtr = nullptr;
        comp_treelets_sbt_record.data.treeletBufferPtr = (device::PKDlet*) treelets_data_[pl_idx];
        comp_treelets_sbt_record.data.radius = particles.GetGlobalRadius();
        //comp_treelets_sbt_record.data.hasGlobalRadius = has_global_radius(particles);
        comp_treelets_sbt_record.data.hasColorData = has_color(particles);
        comp_treelets_sbt_record.data.globalColor =
            glm::vec4(particles.GetGlobalColour()[0] / 255.f, particles.GetGlobalColour()[1] / 255.f,
                particles.GetGlobalColour()[2] / 255.f, particles.GetGlobalColour()[3] / 255.f);
        comp_treelets_sbt_record.data.particleCount = p_count;
        //comp_treelets_sbt_record.data.worldBounds = local_boxes_[pl_idx];

        /*if (!has_global_radius(particles)) {
            treelets_sbt_record.data.radiusBufferPtr = (float*) radius_data_[pl_idx];
        }*/
        //comp_treelets_sbt_record.data.radiusBufferPtr = nullptr;
        if (has_color(particles)) {
            comp_treelets_sbt_record.data.colorBufferPtr = (glm::vec4*) color_data_[pl_idx];
        }
        comp_treelets_sbt_records_.push_back(comp_treelets_sbt_record);

        // occlusion stuff
        SBTRecord<device::QTreeletsGeoData> comp_treelets_sbt_record_occlusion;
        OPTIX_CHECK_ERROR(
            optixSbtRecordPackHeader(comp_treelets_occlusion_module_, &comp_treelets_sbt_record_occlusion));

        comp_treelets_sbt_record_occlusion.data = comp_treelets_sbt_record.data;
        comp_treelets_sbt_records_.push_back(comp_treelets_sbt_record_occlusion);


        SBTRecord<device::STreeletsGeoData> s_comp_treelets_sbt_record;
        OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(s_comp_treelets_module_, &s_comp_treelets_sbt_record));

        s_comp_treelets_sbt_record.data.particleBufferPtr = (device::SPKDParticle*) particle_data_[pl_idx];
        //comp_treelets_sbt_record.data.radiusBufferPtr = nullptr;
        s_comp_treelets_sbt_record.data.colorBufferPtr = nullptr;
        s_comp_treelets_sbt_record.data.treeletBufferPtr = (device::SPKDlet*) treelets_data_[pl_idx];
        s_comp_treelets_sbt_record.data.radius = particles.GetGlobalRadius();
        //comp_treelets_sbt_record.data.hasGlobalRadius = has_global_radius(particles);
        s_comp_treelets_sbt_record.data.hasColorData = has_color(particles);
        s_comp_treelets_sbt_record.data.globalColor =
            glm::vec4(particles.GetGlobalColour()[0] / 255.f, particles.GetGlobalColour()[1] / 255.f,
                particles.GetGlobalColour()[2] / 255.f, particles.GetGlobalColour()[3] / 255.f);
        s_comp_treelets_sbt_record.data.particleCount = p_count;
        //comp_treelets_sbt_record.data.worldBounds = local_boxes_[pl_idx];

        /*if (!has_global_radius(particles)) {
            treelets_sbt_record.data.radiusBufferPtr = (float*) radius_data_[pl_idx];
        }*/
        //comp_treelets_sbt_record.data.radiusBufferPtr = nullptr;
        if (has_color(particles)) {
            s_comp_treelets_sbt_record.data.colorBufferPtr = (glm::vec4*) color_data_[pl_idx];
        }
        s_comp_treelets_sbt_records_.push_back(s_comp_treelets_sbt_record);

        // occlusion stuff
        SBTRecord<device::STreeletsGeoData> s_comp_treelets_sbt_record_occlusion;
        OPTIX_CHECK_ERROR(
            optixSbtRecordPackHeader(s_comp_treelets_occlusion_module_, &s_comp_treelets_sbt_record_occlusion));

        s_comp_treelets_sbt_record_occlusion.data = s_comp_treelets_sbt_record.data;
        s_comp_treelets_sbt_records_.push_back(s_comp_treelets_sbt_record_occlusion);


        SBTRecord<device::QPKDTreeletsGeoData> qpkd_treelets_sbt_record;
        OPTIX_CHECK_ERROR(optixSbtRecordPackHeader(qpkd_treelets_module_, &qpkd_treelets_sbt_record));

        qpkd_treelets_sbt_record.data.particleBufferPtr = (void*) particle_data_[pl_idx];
        qpkd_treelets_sbt_record.data.colorBufferPtr = nullptr;
        qpkd_treelets_sbt_record.data.treeletBufferPtr = (device::QPKDlet*) treelets_data_[pl_idx];
        qpkd_treelets_sbt_record.data.radius = particles.GetGlobalRadius();
        qpkd_treelets_sbt_record.data.hasColorData = has_color(particles);
        qpkd_treelets_sbt_record.data.globalColor =
            glm::vec4(particles.GetGlobalColour()[0] / 255.f, particles.GetGlobalColour()[1] / 255.f,
                particles.GetGlobalColour()[2] / 255.f, particles.GetGlobalColour()[3] / 255.f);
        qpkd_treelets_sbt_record.data.particleCount = p_count;
        qpkd_treelets_sbt_record.data.expXBuffer = (char*) exp_x_data_[pl_idx];
        qpkd_treelets_sbt_record.data.expYBuffer = (char*) exp_y_data_[pl_idx];
        qpkd_treelets_sbt_record.data.expZBuffer = (char*) exp_z_data_[pl_idx];
        qpkd_treelets_sbt_record.data.selectedType = qtreelet_type_slot_.Param<core::param::EnumParam>()->Value();

        if (has_color(particles)) {
            qpkd_treelets_sbt_record.data.colorBufferPtr = (glm::vec4*) color_data_[pl_idx];
        }
        qpkd_treelets_sbt_records_.push_back(qpkd_treelets_sbt_record);

        // occlusion stuff
        SBTRecord<device::QPKDTreeletsGeoData> qpkd_treelets_sbt_record_occlusion;
        OPTIX_CHECK_ERROR(
            optixSbtRecordPackHeader(qpkd_treelets_occlusion_module_, &qpkd_treelets_sbt_record_occlusion));

        qpkd_treelets_sbt_record_occlusion.data = qpkd_treelets_sbt_record.data;
        qpkd_treelets_sbt_records_.push_back(qpkd_treelets_sbt_record_occlusion);
    }

    ++sbt_version;

    return true;
}

} // namespace megamol::optix_hpg
