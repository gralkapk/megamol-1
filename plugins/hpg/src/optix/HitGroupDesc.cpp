#include "stdafx.h"
#include "HitGroupDesc.h"


megamol::hpg::optix::HitGroupDesc::HitGroupDesc()
    : _out_desc_slot("outDesc", ""), _in_shading_slot("inShading", ""), _in_intersection_slot("inIntersection", "") {
    _out_desc_slot.SetCallback(
        CallProgramGroupDesc::ClassName(), CallProgramGroupDesc::FunctionName(0), &HitGroupDesc::get_desc_cb);
    MakeSlotAvailable(&_out_desc_slot);

    _in_shading_slot.SetCompatibleCall<CallProgramGroupDescDescription>();
    MakeSlotAvailable(&_in_shading_slot);

    _in_intersection_slot.SetCompatibleCall<CallProgramGroupDescDescription>();
    MakeSlotAvailable(&_in_intersection_slot);
}


megamol::hpg::optix::HitGroupDesc::~HitGroupDesc() { this->Release(); }


bool megamol::hpg::optix::HitGroupDesc::create() { return true; }


void megamol::hpg::optix::HitGroupDesc::release() {}


bool megamol::hpg::optix::HitGroupDesc::get_desc_cb(core::Call& c) {
    auto out_call = dynamic_cast<CallProgramGroupDesc*>(&c);
    if (out_call != nullptr) return false;

    auto in_shading_call = _in_shading_slot.CallAs<CallProgramGroupDesc>();
    if (in_shading_call != nullptr) return false;

    auto in_intersection_call = _in_intersection_slot.CallAs<CallProgramGroupDesc>();
    if (in_intersection_call != nullptr) return false;

    if (in_shading_call->is_dirty() || in_intersection_call->is_dirty()) {
        auto shading_group = in_shading_call->get_descriptor();
        if (shading_group == nullptr) return false;
        auto intersection_group = in_intersection_call->get_descriptor();
        if (intersection_group == nullptr) return false;

        _desc = std::make_shared<OptixProgramGroupDesc>();
        _desc->hitgroup.moduleCH = shading_group->hitgroup.moduleCH;
        _desc->hitgroup.entryFunctionNameCH = shading_group->hitgroup.entryFunctionNameCH;
        _desc->hitgroup.moduleAH = shading_group->hitgroup.moduleAH;
        _desc->hitgroup.entryFunctionNameAH = shading_group->hitgroup.entryFunctionNameAH;
        _desc->hitgroup.moduleIS = intersection_group->hitgroup.moduleIS;
        _desc->hitgroup.entryFunctionNameIS = intersection_group->hitgroup.entryFunctionNameIS;

        in_shading_call->reset_dirty();
        in_intersection_call->reset_dirty();

        out_call->set_dirty();
    }

    out_call->set_descriptor(_desc);

    return true;
}
