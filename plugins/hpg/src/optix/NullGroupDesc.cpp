#include "stdafx.h"
#include "NullGroupDesc.h"


megamol::hpg::optix::NullGroupDesc::NullGroupDesc() : _out_desc_slot("outDesc", "") {
    _out_desc_slot.SetCallback(
        CallProgramGroupDesc::ClassName(), CallProgramGroupDesc::FunctionName(0), &NullGroupDesc::get_desc_cb);
    MakeSlotAvailable(&_out_desc_slot);

    _desc = std::make_shared<OptixProgramGroupDesc>(0);
}


megamol::hpg::optix::NullGroupDesc::~NullGroupDesc() { this->Release(); }


bool megamol::hpg::optix::NullGroupDesc::create() { return true; }


void megamol::hpg::optix::NullGroupDesc::release() {}


bool megamol::hpg::optix::NullGroupDesc::get_desc_cb(core::Call& c) {
    auto out_call = dynamic_cast<CallProgramGroupDesc*>(&c);
    if (out_call != nullptr) return false;

    out_call->set_descriptor(_desc);

    return true;
}
