#include "Accumulator.h"

#include "mmcore/param/IntParam.h"
#include "mmcore/param/EnumParam.h"

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

        std::vector<std::vector<glm::vec3>> avg_pos;
        std::vector<std::vector<glm::vec3>> avg_dir;
        std::vector<std::vector<glm::vec4>> avg_col;
        std::vector<std::vector<uint64_t>> base_id;

        switch (direction) {
        case dir_t::BACKWARD:
            std::tie(base_id, avg_pos, avg_dir, avg_col) = collect_backward(in_data, window_size);
            break;
        case dir_t::CENTRAL:
            collect_central(in_data, window_size);
            break;
        case dir_t::FORWARD:
            collect_forward(in_data, window_size);
            break;
        }

        frame_id_ = in_data->FrameID();
        in_data_hash_ = in_data->DataHash();
        resetDirty();
	}


    return true;
}


bool megamol::moldyn::Accumulator::get_extent_cb(core::Call& c) {
    return true;
}


std::tuple<std::vector<std::vector<uint64_t>>, std::vector<std::vector<glm::vec3>>, std::vector<std::vector<glm::vec3>>,
    std::vector<std::vector<glm::vec4>>>
megamol::moldyn::Accumulator::collect_backward(geocalls::MultiParticleDataCall* data, int window_size) {
    auto const frame_count = data->FrameCount();
    auto const base_frame_id = static_cast<int>(data->FrameID());

    auto const a_f_id = base_frame_id - window_size < 0 ? 0 : base_frame_id - window_size;
    auto const b_f_id = base_frame_id;

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

    for (int f_id = b_f_id; f_id >= a_f_id; --f_id) {
        do {
            data->SetFrameID(f_id);
            (*data)(0);
        } while (data->FrameID() != f_id);

        auto const pl_count = data->GetParticleListCount();

        if (f_id == base_frame_id) {
            //id_map.resize(pl_count);
            for (std::decay_t<decltype(pl_count)> pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
                auto const& particles = data->AccessParticles(pl_idx);
                auto const p_count = particles.GetCount();
                auto& s_id_map = id_map[pl_idx];
                s_id_map.reserve(p_count);
                avg_pos[pl_idx].resize(p_count, glm::vec3(0));
                avg_dir[pl_idx].resize(p_count, glm::vec3(0));
                avg_col[pl_idx].resize(p_count, glm::vec4(0));
                base_id[pl_idx].resize(p_count);
                auto& s_base_id = base_id[pl_idx];
                auto const id_acc = particles.GetParticleStore().GetIDAcc();
                for (std::decay_t<decltype(p_count)> p_idx = 0; p_idx < p_count; ++p_idx) {
                    auto const id = id_acc->Get_u64(p_idx);
                    s_id_map[id] = p_idx;
                    s_base_id[p_idx] = id;
                }
            }
        }

        for (std::decay_t<decltype(pl_count)> pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
            auto const& particles = data->AccessParticles(pl_idx);
            auto const p_count = particles.GetCount();
            auto const& s_id_map = id_map[pl_idx];
            auto& s_avg_pos = avg_pos[pl_idx];
            auto& s_avg_dir = avg_dir[pl_idx];
            auto& s_avg_col = avg_col[pl_idx];

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
                    s_avg_col[idx] += glm::vec4(cr_acc->Get_f(p_idx), cg_acc->Get_f(p_idx), cb_acc->Get_f(p_idx), ca_acc->Get_f(p_idx));
                }
            }
        }
    }

    for (std::decay_t<decltype(pl_count)> pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
        std::transform(avg_pos[pl_idx].begin(), avg_pos[pl_idx].end(), avg_pos[pl_idx].begin(),
            [&window_size](auto const& val) { return val / static_cast<float>(window_size); });
        std::transform(avg_dir[pl_idx].begin(), avg_dir[pl_idx].end(), avg_dir[pl_idx].begin(),
            [&window_size](auto const& val) { return val / static_cast<float>(window_size); });
        std::transform(avg_col[pl_idx].begin(), avg_col[pl_idx].end(), avg_col[pl_idx].begin(),
            [&window_size](auto const& val) { return val / static_cast<float>(window_size); });
    }

    return std::make_tuple(base_id, avg_pos, avg_dir, avg_col);
}


void megamol::moldyn::Accumulator::collect_central(geocalls::MultiParticleDataCall* data, int window_size) {}


void megamol::moldyn::Accumulator::collect_forward(geocalls::MultiParticleDataCall* data, int window_size) {}
