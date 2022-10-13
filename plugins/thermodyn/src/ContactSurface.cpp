#include "ContactSurface.h"

#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"

#include "geometry_calls/MultiParticleDataCall.h"

#include <glm/glm.hpp>


megamol::thermodyn::ContactSurface::ContactSurface()
        : out_data_slot_("dataOut", "")
        , in_data_slot_("dataIn", "")
        , knn_slot_("#neighbors", "")
        , distance_threshold_slot_("distance threshold", "")
        , min_threshold_slot_("threshold::min", "")
        , max_threshold_slot_("threshold::max", "") {
    out_data_slot_.SetCallback(geocalls::MultiParticleDataCall::ClassName(),
        geocalls::MultiParticleDataCall::FunctionName(0), &ContactSurface::get_data_cb);
    out_data_slot_.SetCallback(geocalls::MultiParticleDataCall::ClassName(),
        geocalls::MultiParticleDataCall::FunctionName(1), &ContactSurface::get_extent_cb);
    MakeSlotAvailable(&out_data_slot_);

    in_data_slot_.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&in_data_slot_);

    knn_slot_ << new core::param::IntParam(1, 1);
    MakeSlotAvailable(&knn_slot_);

    distance_threshold_slot_ << new core::param::FloatParam(3.0f, std::numeric_limits<float>::min());
    MakeSlotAvailable(&distance_threshold_slot_);

    min_threshold_slot_ << new core::param::FloatParam(0.1f, std::numeric_limits<float>::min(), 1.0f);
    MakeSlotAvailable(&min_threshold_slot_);
    max_threshold_slot_ << new core::param::FloatParam(0.5f, std::numeric_limits<float>::min(), 1.0f);
    MakeSlotAvailable(&max_threshold_slot_);
}


megamol::thermodyn::ContactSurface::~ContactSurface() {
    this->Release();
}


bool megamol::thermodyn::ContactSurface::create() {
    return true;
}


void megamol::thermodyn::ContactSurface::release() {}


bool megamol::thermodyn::ContactSurface::get_data_cb(core::Call& c) {
    auto out_call = dynamic_cast<geocalls::MultiParticleDataCall*>(&c);
    if (out_call == nullptr)
        return false;

    auto in_call = in_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (in_call == nullptr)
        return false;

    if (!(*in_call)(0))
        return false;

    if (in_data_hash_ != in_call->DataHash() || frame_id_ != in_call->FrameID() || check_dirty()) {
        auto const pl_count = in_call->GetParticleListCount();
        // expecting two list entries
        if (pl_count != 2) {
            core::utility::log::Log::DefaultLog.WriteError("[ContactSurface] Expecting particle list with two entries");
            return false;
        }
        /*for (std::decay_t<decltype(pl_count)> pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
            auto const& part = in_call->AccessParticles(pl_idx);
        }*/

        auto const& phase_0 = in_call->AccessParticles(0);
        auto const& phase_1 = in_call->AccessParticles(1);

        auto const bbox = in_call->AccessBoundingBoxes().ObjectSpaceBBox();

        std::array<float, 6> a_bbox = {bbox.Left(), bbox.Right(), bbox.Bottom(), bbox.Top(), bbox.Back(), bbox.Front()};

        auto const p1_count = phase_1.GetCount();

        auto const& p1_x_acc = phase_1.GetParticleStore().GetXAcc();
        auto const& p1_y_acc = phase_1.GetParticleStore().GetYAcc();
        auto const& p1_z_acc = phase_1.GetParticleStore().GetZAcc();

        if (in_data_hash_ != in_call->DataHash() || frame_id_ != in_call->FrameID()) {
            std::vector<float> phase_1_pos(p1_count * 3);
            for (std::decay_t<decltype(p1_count)> p_idx = 0; p_idx < p1_count; ++p_idx) {
                phase_1_pos[p_idx * 3 + 0] = p1_x_acc->Get_f(p_idx);
                phase_1_pos[p_idx * 3 + 1] = p1_y_acc->Get_f(p_idx);
                phase_1_pos[p_idx * 3 + 2] = p1_z_acc->Get_f(p_idx);
            }

            std::array<float, 3> weights = {1.f, 1.f, 1.f};

            myPts = std::make_shared<datatools::genericPointcloud<float, 3>>(phase_1_pos, a_bbox, weights);
            core::utility::log::Log::DefaultLog.WriteInfo("[ContactSurface] Creating kd tree");
            particleTree = std::make_shared<my_kd_tree_t>(
                3 /* dim */, *myPts, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
            particleTree->buildIndex();
            core::utility::log::Log::DefaultLog.WriteInfo("[ContactSurface] Created kd tree");
        }

        auto const p0_count = phase_0.GetCount();

        std::list<std::tuple<size_t, size_t, float>> edges;

        auto const num_neighbors = knn_slot_.Param<core::param::IntParam>()->Value();
        auto const distance_threshold = distance_threshold_slot_.Param<core::param::FloatParam>()->Value();
        auto const min_threshold = min_threshold_slot_.Param<core::param::FloatParam>()->Value();
        auto const max_threshold = max_threshold_slot_.Param<core::param::FloatParam>()->Value();

        float query_v[3];
        auto const& p0_x_acc = phase_0.GetParticleStore().GetXAcc();
        auto const& p0_y_acc = phase_0.GetParticleStore().GetYAcc();
        auto const& p0_z_acc = phase_0.GetParticleStore().GetZAcc();
        std::vector<size_t> knn_indices(num_neighbors);
        std::vector<float> knn_distances(num_neighbors);
        for (std::decay_t<decltype(p0_count)> p_idx = 0; p_idx < p0_count; ++p_idx) {
            query_v[0] = p0_x_acc->Get_f(p_idx);
            query_v[1] = p0_y_acc->Get_f(p_idx);
            query_v[2] = p0_z_acc->Get_f(p_idx);
            auto const N = particleTree->knnSearch(query_v, num_neighbors, knn_indices.data(), knn_distances.data());
            for (std::decay_t<decltype(N)> n_idx = 0; n_idx < N; ++n_idx) {
                if (knn_distances[n_idx] <= distance_threshold * 0.25f) {
                    edges.push_back(std::make_tuple(p_idx, knn_indices[n_idx], knn_distances[n_idx]));
                }
            }
        }

        out_data_vec_.clear();
        out_data_vec_.reserve(edges.size());
        out_normal_vec_.clear();
        out_normal_vec_.reserve(edges.size());
        for (auto const& e : edges) {
            auto const s_idx = std::get<0>(e);
            auto const e_idx = std::get<1>(e);
            auto const s_point = glm::vec3(p0_x_acc->Get_f(s_idx), p0_y_acc->Get_f(s_idx), p0_z_acc->Get_f(s_idx));
            auto const e_point = glm::vec3(p1_x_acc->Get_f(e_idx), p1_y_acc->Get_f(e_idx), p1_z_acc->Get_f(e_idx));
            auto const f_point = s_point + 0.5f * (e_point - s_point);
            auto const normal = glm::normalize(e_point - s_point);
            out_data_vec_.emplace_back(f_point.x, f_point.y, f_point.z);
            /*out_data_vec_.push_back(f_point.y);
            out_data_vec_.push_back(f_point.z);*/
            out_normal_vec_.emplace_back(normal.x, normal.y, normal.z);
            /*out_normal_vec_.push_back(normal.y);
            out_normal_vec_.push_back(normal.z);*/
        }

        // clean up data
        // remove high density points
        // remove low density points
        std::vector<size_t> densities(edges.size());
        auto params = nanoflann::SearchParams();
        params.sorted = false;
        for (uint64_t idx = 0; idx < edges.size(); ++idx) {
            query_v[0] = out_data_vec_[idx].x;
            query_v[1] = out_data_vec_[idx].y;
            query_v[2] = out_data_vec_[idx].z;
            std::vector<std::pair<size_t, float>> temp_dis;
            auto const N = particleTree->radiusSearch(query_v, 2.5f, temp_dis, params);
            densities[idx] = N;
        }

        auto minmax_el = std::minmax_element(densities.begin(), densities.end());
        core::utility::log::Log::DefaultLog.WriteInfo(
            "[ContactSurface] Min %d, Max %d", *minmax_el.first, *minmax_el.second);

        auto const range = *minmax_el.second - *minmax_el.first;
        auto const min_thres = *minmax_el.first + min_threshold * range;
        auto const max_thres = *minmax_el.second - max_threshold * range;

        decltype(out_data_vec_) temp_out_data;
        temp_out_data.reserve(out_data_vec_.size());
        decltype(out_normal_vec_) temp_out_normals;
        temp_out_normals.reserve(out_normal_vec_.size());

        for (uint64_t idx = 0; idx < edges.size(); ++idx) {
            if (min_thres < densities[idx] && densities[idx] < max_thres) {
                temp_out_data.emplace_back(out_data_vec_[idx]);
                temp_out_normals.emplace_back(out_normal_vec_[idx]);
            }
        }

        out_data_vec_ = temp_out_data;
        out_normal_vec_ = temp_out_normals;

        in_data_hash_ = in_call->DataHash();
        frame_id_ = in_call->FrameID();
        reset_dirty();
        ++out_data_hash_;
    }

    out_call->SetParticleListCount(1);
    out_call->AccessBoundingBoxes() = in_call->AccessBoundingBoxes();
    auto& out_part = out_call->AccessParticles(0);
    out_part.SetCount(out_data_vec_.size());
    out_part.SetGlobalColour(255, 0, 0);
    out_part.SetGlobalRadius(0.5f);
    out_part.SetVertexData(geocalls::SimpleSphericalParticles::VERTDATA_FLOAT_XYZ, out_data_vec_.data());
    out_part.SetDirData(geocalls::SimpleSphericalParticles::DIRDATA_FLOAT_XYZ, out_normal_vec_.data());
    out_call->SetDataHash(out_data_hash_);
    out_call->SetFrameID(frame_id_);

    return true;
}


bool megamol::thermodyn::ContactSurface::get_extent_cb(core::Call& c) {
    auto out_call = dynamic_cast<geocalls::MultiParticleDataCall*>(&c);
    if (out_call == nullptr)
        return false;

    auto in_call = in_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (in_call == nullptr)
        return false;

    in_call->SetFrameID(out_call->FrameID());
    if (!(*in_call)(1))
        return false;

    out_call->AccessBoundingBoxes() = in_call->AccessBoundingBoxes();
    out_call->SetFrameCount(in_call->FrameCount());

    return true;
}
