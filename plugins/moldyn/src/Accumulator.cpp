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

        switch (direction) {
        case dir_t::BACKWARD:
            collect_backward(in_data, window_size);
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


void megamol::moldyn::Accumulator::collect_backward(geocalls::MultiParticleDataCall* data, int window_size) {
    auto const frame_count = data->FrameCount();
    auto const base_frame_id = static_cast<int>(data->FrameID());

    auto const a_f_id = base_frame_id - window_size < 0 ? 0 : base_frame_id - window_size;
    auto const b_f_id = base_frame_id;

    std::vector<std::unordered_map<uint64_t /*id*/, uint64_t /*idx*/>> id_map;

    for (int f_id = b_f_id; f_id >= a_f_id; --f_id) {
        do {
            data->SetFrameID(f_id);
            (*data)(0);
        } while (data->FrameID() != f_id);

        auto const pl_count = data->GetParticleListCount();

        if (f_id == base_frame_id) {
            id_map.resize(pl_count);
            for (std::decay_t<decltype(pl_count)> pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
                id_map[pl_idx].reserve(data->AccessParticles(pl_idx).GetCount());
            }
        }

        for (std::decay_t<decltype(pl_count)> pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
            auto const& particles = data->AccessParticles(pl_idx);
            auto const p_count = particles.GetCount();
            for (std::decay_t<decltype(p_count)> p_idx = 0; p_idx < p_count; ++p_idx) {
            }
        }
    }
}


void megamol::moldyn::Accumulator::collect_central(geocalls::MultiParticleDataCall* data, int window_size) {}


void megamol::moldyn::Accumulator::collect_forward(geocalls::MultiParticleDataCall* data, int window_size) {}
