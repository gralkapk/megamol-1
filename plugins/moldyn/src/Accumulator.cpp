#include "Accumulator.h"

#include "mmcore/param/EnumParam.h"
#include "mmcore/param/IntParam.h"

#include "geometry_calls/MultiParticleDataCall.h"

#include <unordered_map>


megamol::moldyn::Accumulator::Accumulator()
        : data_in_slot_("dataIn", "")
        , data_out_slot_("dataOut", "")
        , window_size_slot_("windowSize", "")
        , direction_slot_("direction", "") {
    data_in_slot_.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&data_in_slot_);

    data_out_slot_.SetCallback(geocalls::MultiParticleDataCall::ClassName(),
        geocalls::MultiParticleDataCall::FunctionName(0), &Accumulator::get_data_cb);
    data_out_slot_.SetCallback(geocalls::MultiParticleDataCall::ClassName(),
        geocalls::MultiParticleDataCall::FunctionName(1), &Accumulator::get_extent_cb);
    MakeSlotAvailable(&data_out_slot_);

    window_size_slot_ << new core::param::IntParam(1, 1);
    MakeSlotAvailable(&window_size_slot_);

    auto ep = new core::param::EnumParam(static_cast<int>(dir_t::BACKWARD));
    ep->SetTypePair(static_cast<int>(dir_t::BACKWARD), "Backward");
    ep->SetTypePair(static_cast<int>(dir_t::CENTRAL), "Central");
    ep->SetTypePair(static_cast<int>(dir_t::FORWARD), "Forward");
    direction_slot_ << ep;
    MakeSlotAvailable(&direction_slot_);
}


megamol::moldyn::Accumulator::~Accumulator() {
    this->Release();
}


bool megamol::moldyn::Accumulator::create() {
    return true;
}


void megamol::moldyn::Accumulator::release() {}


bool megamol::moldyn::Accumulator::get_data_cb(core::Call& c) {
    auto out_data = dynamic_cast<geocalls::MultiParticleDataCall*>(&c);
    if (out_data == nullptr)
        return false;

    auto in_data = data_in_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (in_data == nullptr)
        return false;

    in_data->SetFrameID(out_data->FrameID());
    if (!(*in_data)(0))
        return false;

    if (in_data->FrameID() != frame_id_ || in_data->DataHash() != in_data_hash_ || isDirty()) {
        auto const window_size = window_size_slot_.Param<core::param::IntParam>()->Value();
        auto const direction = static_cast<dir_t>(direction_slot_.Param<core::param::EnumParam>()->Value());

        auto const tmp_fid = in_data->FrameID();
        auto const tmp_dh = in_data->DataHash();

        auto const pl_count = in_data->GetParticleListCount();
        global_radii_.resize(pl_count);
        base_pos_.resize(pl_count);
        for (std::decay_t<decltype(pl_count)> pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
            auto const& particles = in_data->AccessParticles(pl_idx);
            global_radii_[pl_idx] = particles.GetGlobalRadius();
            auto const p_count = particles.GetCount();
            auto& s_pos = base_pos_[pl_idx];
            s_pos.resize(particles.GetCount());
            auto const x_acc = particles.GetParticleStore().GetXAcc();
            auto const y_acc = particles.GetParticleStore().GetYAcc();
            auto const z_acc = particles.GetParticleStore().GetZAcc();
            for (std::decay_t<decltype(p_count)> p_idx = 0; p_idx < p_count; ++p_idx) {
                s_pos[p_idx] = glm::vec3(x_acc->Get_f(p_idx), y_acc->Get_f(p_idx), z_acc->Get_f(p_idx));
            }
        }

        std::vector<std::vector<glm::vec3>> avg_pos;
        std::vector<std::vector<glm::vec3>> avg_dir;
        std::vector<std::vector<glm::vec4>> avg_col;
        std::vector<std::vector<uint64_t>> base_id;

        switch (direction) {
        case dir_t::BACKWARD:
            std::tie(base_id, avg_pos, avg_dir, avg_col) = collect_backward(in_data, window_size);
            break;
        case dir_t::CENTRAL:
            std::tie(base_id, avg_pos, avg_dir, avg_col) = collect_central(in_data, window_size);
            break;
        case dir_t::FORWARD:
            std::tie(base_id, avg_pos, avg_dir, avg_col) = collect_forward(in_data, window_size);
            break;
        }

        avg_pos_ = avg_pos;
        avg_dir_ = avg_dir;
        avg_col_ = avg_col;
        base_id_ = base_id;

        frame_id_ = tmp_fid;
        in_data_hash_ = tmp_dh;
        resetDirty();
        ++out_data_hash_;
    }

    out_data->SetParticleListCount(avg_pos_.size());
    for (unsigned int pl_idx = 0; pl_idx < avg_pos_.size(); ++pl_idx) {
        auto& particles = out_data->AccessParticles(pl_idx);
        auto const& s_pos = base_pos_[pl_idx];
        auto const& s_dir = avg_dir_[pl_idx];
        auto const& s_col = avg_col_[pl_idx];
        auto const& s_base_id = base_id_[pl_idx];
        particles.SetCount(s_pos.size());
        if (s_pos.size() != 0) {
            particles.SetVertexData(geocalls::SimpleSphericalParticles::VERTDATA_FLOAT_XYZ, s_pos.data());
            particles.SetDirData(geocalls::SimpleSphericalParticles::DIRDATA_FLOAT_XYZ, s_dir.data());
            particles.SetColourData(geocalls::SimpleSphericalParticles::COLDATA_FLOAT_RGBA, s_col.data());
            particles.SetIDData(geocalls::SimpleSphericalParticles::IDDATA_UINT64, s_base_id.data());
        }
        particles.SetGlobalRadius(global_radii_[pl_idx]);
    }
    out_data->SetFrameID(frame_id_);
    out_data->SetDataHash(out_data_hash_);

    return true;
}


bool megamol::moldyn::Accumulator::get_extent_cb(core::Call& c) {
    auto out_data = dynamic_cast<geocalls::MultiParticleDataCall*>(&c);
    if (out_data == nullptr)
        return false;

    auto in_data = data_in_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (in_data == nullptr)
        return false;

    in_data->SetFrameID(out_data->FrameID());
    if (!(*in_data)(1))
        return false;

    out_data->SetFrameCount(in_data->FrameCount());
    out_data->AccessBoundingBoxes() = in_data->AccessBoundingBoxes();

    return true;
}


std::tuple<std::vector<std::vector<uint64_t>>, std::vector<std::vector<glm::vec3>>, std::vector<std::vector<glm::vec3>>,
    std::vector<std::vector<glm::vec4>>>
megamol::moldyn::Accumulator::compute_collection(
    geocalls::MultiParticleDataCall* data, int base_frame_id, int a_f_id, int b_f_id) {
    auto const ws = b_f_id - a_f_id;

    auto const pl_count = data->GetParticleListCount();
    std::vector<std::unordered_map<uint64_t /*id*/, uint64_t /*idx*/>> id_map;
    id_map.resize(pl_count);
    std::vector<std::vector<glm::vec3>> avg_pos;
    avg_pos.resize(pl_count);
    std::vector<std::vector<glm::vec3>> avg_dir;
    avg_dir.resize(pl_count);
    std::vector<std::vector<glm::vec4>> avg_col;
    avg_col.resize(pl_count);
    std::vector<std::vector<uint64_t>> base_id;
    base_id.resize(pl_count);
    std::vector<std::vector<float>> div;
    div.resize(pl_count);

    for (std::decay_t<decltype(pl_count)> pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
        auto const& particles = data->AccessParticles(pl_idx);
        auto const p_count = particles.GetCount();
        auto& s_id_map = id_map[pl_idx];
        s_id_map.clear();
        s_id_map.reserve(p_count);
        avg_pos[pl_idx].resize(p_count, glm::vec3(0));
        avg_dir[pl_idx].resize(p_count, glm::vec3(0));
        avg_col[pl_idx].resize(p_count, glm::vec4(0));
        base_id[pl_idx].resize(p_count);
        div[pl_idx].resize(p_count, 0);
        auto& s_base_id = base_id[pl_idx];
        auto const id_acc = particles.GetParticleStore().GetIDAcc();
        for (std::decay_t<decltype(p_count)> p_idx = 0; p_idx < p_count; ++p_idx) {
            auto const id = id_acc->Get_u64(p_idx);
            s_id_map[id] = p_idx;
            s_base_id[p_idx] = id;
        }
    }

    for (int f_id = b_f_id; f_id >= a_f_id; --f_id) {
        do {
            data->SetFrameID(f_id);
            (*data)(0);
        } while (data->FrameID() != f_id);

        for (std::decay_t<decltype(pl_count)> pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
            auto const& particles = data->AccessParticles(pl_idx);
            auto const p_count = particles.GetCount();
            auto const& s_id_map = id_map[pl_idx];
            auto& s_avg_pos = avg_pos[pl_idx];
            auto& s_avg_dir = avg_dir[pl_idx];
            auto& s_avg_col = avg_col[pl_idx];
            auto& s_div = div[pl_idx];

            auto const id_acc = particles.GetParticleStore().GetIDAcc();

            auto const x_acc = particles.GetParticleStore().GetXAcc();
            auto const y_acc = particles.GetParticleStore().GetYAcc();
            auto const z_acc = particles.GetParticleStore().GetZAcc();

            auto const dx_acc = particles.GetParticleStore().GetDXAcc();
            auto const dy_acc = particles.GetParticleStore().GetDYAcc();
            auto const dz_acc = particles.GetParticleStore().GetDZAcc();

            auto const cr_acc = particles.GetParticleStore().GetCRAcc();
            auto const cg_acc = particles.GetParticleStore().GetCGAcc();
            auto const cb_acc = particles.GetParticleStore().GetCBAcc();
            auto const ca_acc = particles.GetParticleStore().GetCAAcc();

            for (std::decay_t<decltype(p_count)> p_idx = 0; p_idx < p_count; ++p_idx) {
                auto const id = id_acc->Get_u64(p_idx);
                auto const fit = s_id_map.find(id);
                if (fit != s_id_map.end()) {
                    // we have found a match
                    auto const idx = fit->second;

                    s_avg_pos[idx] += glm::vec3(x_acc->Get_f(p_idx), y_acc->Get_f(p_idx), z_acc->Get_f(p_idx));
                    s_avg_dir[idx] += glm::vec3(dx_acc->Get_f(p_idx), dy_acc->Get_f(p_idx), dz_acc->Get_f(p_idx));
                    s_avg_col[idx] += glm::vec4(
                        cr_acc->Get_f(p_idx), cg_acc->Get_f(p_idx), cb_acc->Get_f(p_idx), ca_acc->Get_f(p_idx));
                    s_div[idx] += 1.f;
                }
            }
        }
    }

    for (std::decay_t<decltype(pl_count)> pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
        std::transform(avg_pos[pl_idx].begin(), avg_pos[pl_idx].end(), div[pl_idx].begin(), avg_pos[pl_idx].begin(),
            [&ws](auto const& val, auto const& d) { return val / d; });
        std::transform(avg_dir[pl_idx].begin(), avg_dir[pl_idx].end(), div[pl_idx].begin(), avg_dir[pl_idx].begin(),
            [&ws](auto const& val, auto const& d) { return val / d; });
        std::transform(avg_col[pl_idx].begin(), avg_col[pl_idx].end(), div[pl_idx].begin(), avg_col[pl_idx].begin(),
            [&ws](auto const& val, auto const& d) { return val / d; });
    }

    return std::make_tuple(base_id, avg_pos, avg_dir, avg_col);
}


std::tuple<std::vector<std::vector<uint64_t>>, std::vector<std::vector<glm::vec3>>, std::vector<std::vector<glm::vec3>>,
    std::vector<std::vector<glm::vec4>>>
megamol::moldyn::Accumulator::collect_backward(geocalls::MultiParticleDataCall* data, int window_size) {
    auto const frame_count = data->FrameCount();
    auto const base_frame_id = static_cast<int>(data->FrameID());

    auto const a_f_id = base_frame_id - window_size < 0 ? 0 : base_frame_id - window_size;
    auto const b_f_id = base_frame_id;

    return compute_collection(data, base_frame_id, a_f_id, b_f_id);
}


std::tuple<std::vector<std::vector<uint64_t>>, std::vector<std::vector<glm::vec3>>, std::vector<std::vector<glm::vec3>>,
    std::vector<std::vector<glm::vec4>>>
megamol::moldyn::Accumulator::collect_central(geocalls::MultiParticleDataCall* data, int window_size) {
    auto const frame_count = data->FrameCount();
    auto const base_frame_id = static_cast<int>(data->FrameID());

    auto const low = static_cast<int>(std::floorf(static_cast<float>(window_size) * 0.5f));
    auto const high = window_size - low;

    auto const a_f_id = base_frame_id - low < 0 ? 0 : base_frame_id - low;
    auto const b_f_id = base_frame_id + high >= frame_count ? frame_count - 1 : base_frame_id + high;

    return compute_collection(data, base_frame_id, a_f_id, b_f_id);
}


std::tuple<std::vector<std::vector<uint64_t>>, std::vector<std::vector<glm::vec3>>, std::vector<std::vector<glm::vec3>>,
    std::vector<std::vector<glm::vec4>>>
megamol::moldyn::Accumulator::collect_forward(geocalls::MultiParticleDataCall* data, int window_size) {
    auto const frame_count = data->FrameCount();
    auto const base_frame_id = static_cast<int>(data->FrameID());

    auto const a_f_id = base_frame_id;
    auto const b_f_id = base_frame_id + window_size >= frame_count ? frame_count - 1 : base_frame_id + window_size;

    return compute_collection(data, base_frame_id, a_f_id, b_f_id);
}
