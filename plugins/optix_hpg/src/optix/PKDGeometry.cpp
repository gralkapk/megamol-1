#include "PKDGeometry.h"

#include "mmcore/param/BoolParam.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/IntParam.h"

namespace megamol::optix_hpg {
PKDGeometry::PKDGeometry()
        : out_geo_slot_("outGeo", "")
        , in_data_slot_("inData", "")
        , mode_slot_("mode", "")
        , compression_slot_("compression", "")
        , threshold_slot_("threshold", "") {
    in_data_slot_.SetCompatibleCall<geocalls::MultiParticleDataCall>();
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

void PKDGeometry::release() {}

bool PKDGeometry::get_data_cb(core::Call& c) {
    return true;
}

bool PKDGeometry::get_extents_cb(core::Call& c) {
    return true;
}

bool PKDGeometry::init(Context const& ctx) {
    return true;
}

bool PKDGeometry::assert_data(geocalls::MultiParticleDataCall const& call, Context const& ctx) {
    free_data(ctx);

    auto const pl_count = call.GetParticleListCount();



    return true;
}
} // namespace megamol::optix_hpg
