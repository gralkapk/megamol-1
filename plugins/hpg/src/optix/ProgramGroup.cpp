#include "stdafx.h"
#include "ProgramGroup.h"

#include "hpg/optix/CallContext.h"
#include "hpg/optix/CallProgramGroup.h"


megamol::hpg::optix::ProgramGroup::ProgramGroup()
    : _out_groups_slot("outGroups", "")
    , _in_context_slot("inContext", "")
    , _in_raygen_slot("inRaygen", "")
    , _in_miss_slot("inMiss", "")
    , _in_exception_slot("inException", "")
    , _in_hitgroup_slot("inHitgroup", "")
    , _in_callables_slot("inCallables", "") {
    _out_groups_slot.SetCallback(
        CallProgramGroupDesc::ClassName(), CallProgramGroupDesc::FunctionName(0), &ProgramGroup::get_desc_cb);
    MakeSlotAvailable(&_out_groups_slot);

    _in_context_slot.SetCompatibleCall<CallContextDescription>();
    MakeSlotAvailable(&_in_context_slot);

    _in_raygen_slot.SetCompatibleCall<CallProgramGroupDescDescription>();
    MakeSlotAvailable(&_in_raygen_slot);

    _in_miss_slot.SetCompatibleCall<CallProgramGroupDescDescription>();
    MakeSlotAvailable(&_in_miss_slot);

    _in_exception_slot.SetCompatibleCall<CallProgramGroupDescDescription>();
    MakeSlotAvailable(&_in_exception_slot);

    _in_hitgroup_slot.SetCompatibleCall<CallProgramGroupDescDescription>();
    MakeSlotAvailable(&_in_hitgroup_slot);

    _in_callables_slot.SetCompatibleCall<CallProgramGroupDescDescription>();
    MakeSlotAvailable(&_in_callables_slot);

    _groups = std::make_shared<std::vector<OptixProgramGroup>>();
}


megamol::hpg::optix::ProgramGroup::~ProgramGroup() { this->Release(); }


bool megamol::hpg::optix::ProgramGroup::create() { return true; }


void megamol::hpg::optix::ProgramGroup::release() {}


bool megamol::hpg::optix::ProgramGroup::get_desc_cb(core::Call& c) {
    auto out_call = dynamic_cast<CallProgramGroup*>(&c);
    if (out_call != nullptr) return false;

    auto in_context_call = _in_raygen_slot.CallAs<CallContext>();
    if (in_context_call != nullptr) return false;

    auto in_raygen_call = _in_raygen_slot.CallAs<CallProgramGroupDesc>();
    if (in_raygen_call != nullptr) return false;

    auto in_miss_call = _in_miss_slot.CallAs<CallProgramGroupDesc>();
    if (in_miss_call != nullptr) return false;

    auto in_exception_call = _in_exception_slot.CallAs<CallProgramGroupDesc>();
    if (in_exception_call != nullptr) return false;

    auto in_hitgroup_call = _in_hitgroup_slot.CallAs<CallProgramGroupDesc>();
    if (in_hitgroup_call != nullptr) return false;

    auto in_callables_call = _in_callables_slot.CallAs<CallProgramGroupDesc>();
    if (in_callables_call != nullptr) return false;

    if (in_context_call->is_dirty() || in_raygen_call->is_dirty() || in_miss_call->is_dirty() ||
        in_exception_call->is_dirty() || in_hitgroup_call->is_dirty() || in_callables_call->is_dirty()) {
        auto context = in_context_call->get_context();
        if (context == nullptr) return false;
        auto raygen_group = in_raygen_call->get_descriptor();
        auto miss_group = in_miss_call->get_descriptor();
        auto exception_group = in_exception_call->get_descriptor();
        auto hit_group = in_hitgroup_call->get_descriptor();
        auto callables_group = in_callables_call->get_descriptor();

        std::vector<OptixProgramGroupDesc> desc;
        desc.reserve(5);

        if (raygen_group != nullptr) {
            desc.push_back(*raygen_group);
        }
        if (miss_group != nullptr) {
            desc.push_back(*miss_group);
        }
        if (exception_group != nullptr) {
            desc.push_back(*exception_group);
        }
        if (hit_group != nullptr) {
            desc.push_back(*hit_group);
        }
        if (callables_group != nullptr) {
            desc.push_back(*callables_group);
        }

        _groups->resize(desc.size());

        OptixProgramGroupOptions options = {0};

        optixProgramGroupCreate(*context, desc.data(), desc.size(), &options, nullptr, nullptr, _groups->data());

        in_raygen_call->reset_dirty();
        in_miss_call->reset_dirty();
        in_exception_call->reset_dirty();
        in_hitgroup_call->reset_dirty();
        in_callables_call->reset_dirty();

        in_context_call->reset_dirty();

        out_call->set_dirty();
    }

    out_call->set_groups(_groups);

    return true;
}
