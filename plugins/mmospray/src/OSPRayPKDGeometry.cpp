/*
 * OSPRayPKDGeometry.cpp
 * Copyright (C) 2009-2017 by MegaMol Team
 * Alle Rechte vorbehalten.
 */

#include <type_traits>

#include "OSPRayPKDGeometry.h"
#include "geometry_calls/MultiParticleDataCall.h"
#include "mmcore/Call.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"
#include "mmcore/param/Vector3fParam.h"
#include "mmcore/utility/log/Log.h"
#include "vislib/forceinline.h"

#include "mmospray/CallOSPRayAPIObject.h"
#include "mmstd/renderer/CallClipPlane.h"
#include "ospray/ospray_cpp.h"

#include "rkcommon/math/box.h"
#include "rkcommon/math/vec.h"

#include <tbb/parallel_for.h>

namespace megamol::ospray {


OSPRayPKDGeometry::OSPRayPKDGeometry()
        : getDataSlot("getdata", "Connects to the data source")
        , deployStructureSlot("deployStructureSlot", "Connects to an OSPRayAPIStructure")
        , mode_slot_("mode", "")
        , threshold_slot_("treelets::threshold", "")
/*, colorTypeSlot("colorType", "Set the type of encoded color")*/ {

    this->getDataSlot.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    this->MakeSlotAvailable(&this->getDataSlot);

    auto ep = new core::param::EnumParam(static_cast<int>(mode::PKD));
    ep->SetTypePair(static_cast<int>(mode::PKD), "PKD");
    ep->SetTypePair(static_cast<int>(mode::TREELETS), "TREELETS");
    mode_slot_ << ep;
    MakeSlotAvailable(&mode_slot_);

    threshold_slot_ << new core::param::IntParam(256, 16, 2048);
    MakeSlotAvailable(&threshold_slot_);

    /*auto ep = new megamol::core::param::EnumParam(0);
    ep->SetTypePair(0, "none");
    ep->SetTypePair(1, "RGBu8");
    ep->SetTypePair(2, "RGBAu8");
    ep->SetTypePair(3, "RGBf");
    ep->SetTypePair(4, "RGBAf");
    ep->SetTypePair(5, "I");
    this->colorTypeSlot << ep;
    this->MakeSlotAvailable(&this->colorTypeSlot);*/

    this->deployStructureSlot.SetCallback(
        CallOSPRayAPIObject::ClassName(), CallOSPRayAPIObject::FunctionName(0), &OSPRayPKDGeometry::getDataCallback);
    this->deployStructureSlot.SetCallback(
        CallOSPRayAPIObject::ClassName(), CallOSPRayAPIObject::FunctionName(1), &OSPRayPKDGeometry::getExtendsCallback);
    this->deployStructureSlot.SetCallback(
        CallOSPRayAPIObject::ClassName(), CallOSPRayAPIObject::FunctionName(2), &OSPRayPKDGeometry::getDirtyCallback);
    this->MakeSlotAvailable(&this->deployStructureSlot);
}


bool has_global_color(geocalls::SimpleSphericalParticles::ColourDataType const& ctype) {
    return ctype == geocalls::SimpleSphericalParticles::ColourDataType::COLDATA_NONE ||
           ctype == geocalls::SimpleSphericalParticles::ColourDataType::COLDATA_FLOAT_I ||
           ctype == geocalls::SimpleSphericalParticles::ColourDataType::COLDATA_DOUBLE_I;
}

::ospray::cpp::Geometry set_geometry(geocalls::MultiParticleDataCall::Particles const& parts,
    std::vector<glm::vec3> const& position, std::vector<glm::u8vec4> const& color, datatools::box3f const& bounds, std::vector<datatools::pkdlet> const& treelets) {
    ::ospray::cpp::Geometry geo = ospNewGeometry("PKDGeometry");

    auto const size = position.size();

    auto positionData = ::ospray::cpp::SharedData(position.data(), OSP_VEC3F, size);
    positionData.commit();

    // set bbox
    float box[] = {bounds.lower.x, bounds.lower.y, bounds.lower.z, bounds.upper.x, bounds.upper.y, bounds.upper.z};
    /*rkcommon::math::box3f box;
    box.lower = {bounds.lower.x, bounds.lower.y, bounds.lower.z};
    box.upper = {bounds.upper.x, bounds.upper.y, bounds.upper.z};*/
    auto boundsData = ::ospray::cpp::CopiedData(box, OSP_FLOAT, 6);
    boundsData.commit();

    // set global color
    auto globalColorData = ::ospray::cpp::CopiedData(parts.GetGlobalColour(), OSP_UCHAR, 4);
    globalColorData.commit();

    /* Interface:
    global_radius = getParam<float>("global_radius", 0.5f);

    has_global_color = getParam<bool>("has_global_color", true);
    global_color = getParam<vec4uc>("global_color", vec4uc(255, 0, 0, 255));

    positionData = getParamDataT<vec3f>("position");
    colorData = getParamDataT<vec4uc>("color");

    num_particles = getParam<unsigned int>("num_particles");

    bounds = getParam<box3f>("bounds");
    */

    geo.setParam("global_radius", parts.GetGlobalRadius());
    geo.setParam("has_global_color", has_global_color(parts.GetColourDataType()));
    geo.setParam("global_color", globalColorData);
    geo.setParam("num_particles", static_cast<unsigned int>(size));
    geo.setParam("bounds", boundsData);

    geo.setParam("position", positionData);
    if (!has_global_color(parts.GetColourDataType())) {
        auto colorData = ::ospray::cpp::SharedData(color.data(), OSP_VEC4UC, size);
        colorData.commit();
        geo.setParam("color", colorData);
    }
    if (!treelets.empty()) {
        auto treeletsData =
            ::ospray::cpp::SharedData(treelets.data(), OSP_CHAR, treelets.size() * sizeof(datatools::pkdlet));
        treeletsData.commit();
        geo.setParam("treelets", treeletsData);
    }

    geo.commit();

    return geo;
}

bool OSPRayPKDGeometry::getDataCallback(megamol::core::Call& call) {

    // read Data, calculate  shape parameters, fill data vectors
    auto os = dynamic_cast<CallOSPRayAPIObject*>(&call);
    megamol::geocalls::MultiParticleDataCall* cd = this->getDataSlot.CallAs<megamol::geocalls::MultiParticleDataCall>();


    //auto const minFrameCount = cd->FrameCount();

    //if (minFrameCount == 0) return false;

    //auto frameTime = 0;

    //if (os->FrameID() >= minFrameCount) {
    //    cd->SetFrameID(minFrameCount - 1, true); // isTimeForced flag set to true
    //    frameTime = minFrameCount - 1;
    //} else {
    //    cd->SetFrameID(os->FrameID(), true); // isTimeForced flag set to true
    //    frameTime = os->FrameID();
    //}
    cd->SetFrameID(os->FrameID(), true);
    if (!(*cd)(1))
        return false;

    if (this->datahash != cd->DataHash() || this->time != cd->FrameID() || this->InterfaceIsDirty()) {
        this->datahash = cd->DataHash();
        this->time = cd->FrameID();
    } else {
        return true;
    }

    if (!(*cd)(0))
        return false;

    size_t listCount = cd->GetParticleListCount();

    position.resize(listCount);
    color.resize(listCount);
    treelets.resize(listCount);

    geo_.clear();
    for (size_t i = 0; i < listCount; ++i) {

        geocalls::MultiParticleDataCall::Particles const& parts = cd->AccessParticles(i);
        if (parts.GetCount() == 0)
            continue;

        //auto colorType = this->colorTypeSlot.Param<megamol::core::param::EnumParam>()->Value();

        // TODO build PKD
        auto data = datatools::collectData(parts);
        datatools::box3f bounds;
        std::tie(position[i], color[i], bounds) = data;

        if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(mode::PKD)) {
            datatools::makePKD(position[i], color[i], bounds);
            
            geo_.emplace_back(set_geometry(parts, position[i], color[i], bounds, std::vector<datatools::pkdlet>()));
        } else if (mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(mode::TREELETS)) {
            treelets[i] = datatools::prePartition_inPlace(position[i],
                threshold_slot_.Param<core::param::IntParam>()->Value(), parts.GetGlobalRadius(), color[i]);

            tbb::parallel_for((size_t) 0, treelets[i].size(), [&](size_t treeletID) {
                datatools::makePKD(position[i], color[i], treelets[i][treeletID].bounds, treelets[i][treeletID].begin,
                    treelets[i][treeletID].end);
            });

            for (auto const& t : treelets) {
                geo_.emplace_back(set_geometry(parts, position[i], color[i], bounds, treelets[i]));
            }
        }

#if 0
        geo_.emplace_back(ospNewGeometry("PKDGeometry"));

        auto positionData = ::ospray::cpp::SharedData(position[i].data(), OSP_VEC3F, position[i].size());
        positionData.commit();

        // set bbox
        float box[] = {bounds.lower.x, bounds.lower.y, bounds.lower.z, bounds.upper.x, bounds.upper.y, bounds.upper.z};
        /*rkcommon::math::box3f box;
        box.lower = {bounds.lower.x, bounds.lower.y, bounds.lower.z};
        box.upper = {bounds.upper.x, bounds.upper.y, bounds.upper.z};*/
        auto boundsData = ::ospray::cpp::CopiedData(box, OSP_FLOAT, 6);
        boundsData.commit();

        // set global color
        auto globalColorData = ::ospray::cpp::CopiedData(parts.GetGlobalColour(), OSP_UCHAR, 4);
        globalColorData.commit();

        /* Interface:
        global_radius = getParam<float>("global_radius", 0.5f);

        has_global_color = getParam<bool>("has_global_color", true);
        global_color = getParam<vec4uc>("global_color", vec4uc(255, 0, 0, 255));

        positionData = getParamDataT<vec3f>("position");
        colorData = getParamDataT<vec4uc>("color");

        num_particles = getParam<unsigned int>("num_particles");

        bounds = getParam<box3f>("bounds");
        */

        geo_.back().setParam("global_radius", parts.GetGlobalRadius());
        geo_.back().setParam("has_global_color", has_global_color(parts.GetColourDataType()));
        geo_.back().setParam("global_color", globalColorData);
        geo_.back().setParam("num_particles", static_cast<unsigned int>(parts.GetCount()));
        geo_.back().setParam("bounds", boundsData);

        geo_.back().setParam("position", positionData);
        if (!has_global_color(parts.GetColourDataType())) {
            auto colorData = ::ospray::cpp::SharedData(color[i].data(), OSP_VEC4UC, color[i].size());
            colorData.commit();
            geo_.back().setParam("color", colorData);
        }

        geo_.back().commit();
#endif

        //geo.back().setParam("radius", parts.GetGlobalRadius());
        ////ospSet1i(geo.back(), "colorType", colorType);
        //geo.back().setParam("colorType", 2);
        //geo.back().setParam("position", vertexData);
        //// ospSetData(geo.back(), "bbox", bboxData);
        //geo.back().setParam("bbox", NULL);
        //geo.back().commit();

        // TODO: implement distributed stuff
        // if (this->rd_type.Param<megamol::core::param::EnumParam>()->Value() == MPI_RAYCAST) {
        //    auto const half_radius = element.globalRadius * 0.5f;

        //    auto const bbox = element.boundingBox->ObjectSpaceBBox().PeekBounds(); //< TODO Not all geometries expose
        //    bbox ospcommon::vec3f lower{bbox[0] - half_radius, bbox[1] - half_radius,
        //        bbox[2] - half_radius}; //< TODO The bbox needs to include complete sphere bound
        //    ospcommon::vec3f upper{bbox[3] + half_radius, bbox[4] + half_radius, bbox[5] + half_radius};
        //    // ghostRegions.emplace_back(lower, upper);
        //    worldBounds.extend({lower, upper}); //< TODO Possible hazard if bbox is not centered
        //    regions.emplace_back(lower, upper);
        //}
    }


    std::vector<void*> geo_transfer(geo_.size());
    for (auto i = 0; i < geo_.size(); i++) {
        geo_transfer[i] = geo_[i].handle();
    }
    os->setStructureType(GEOMETRY);
    os->setAPIObjects(std::move(geo_transfer));

    return true;
}


OSPRayPKDGeometry::~OSPRayPKDGeometry() {
    this->Release();
}

bool OSPRayPKDGeometry::create() {
    auto error = ospLoadModule("pkd");
    if (error != OSPError::OSP_NO_ERROR) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "Unable to load OSPRay module: PKD. Error occured in %s:%d", __FILE__, __LINE__);
    }
    return true;
}

void OSPRayPKDGeometry::release() {}

/*
ospray::OSPRayPKDGeometry::InterfaceIsDirty()
*/
bool OSPRayPKDGeometry::InterfaceIsDirty() {
    if (this->mode_slot_.IsDirty() ||
        ((mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(mode::TREELETS)) &&
            this->threshold_slot_.IsDirty())) {
        this->mode_slot_.ResetDirty();
        threshold_slot_.ResetDirty();
        return true;
    } else {
        return false;
    }
}

bool OSPRayPKDGeometry::InterfaceIsDirtyNoReset() const {
    return this->mode_slot_.IsDirty() ||
           ((mode_slot_.Param<core::param::EnumParam>()->Value() == static_cast<int>(mode::TREELETS)) &&
               this->threshold_slot_.IsDirty());
}


bool OSPRayPKDGeometry::getExtendsCallback(core::Call& call) {
    auto os = dynamic_cast<CallOSPRayAPIObject*>(&call);
    geocalls::MultiParticleDataCall* cd = this->getDataSlot.CallAs<geocalls::MultiParticleDataCall>();

    if (cd == NULL)
        return false;
    cd->SetFrameID(os->FrameID(), true); // isTimeForced flag set to true
    // if (!(*cd)(1)) return false; // table returns flase at first attempt and breaks everything
    (*cd)(1);

    core::BoundingBoxes_2 box;
    box.SetBoundingBox(cd->AccessBoundingBoxes().ObjectSpaceBBox());

    os->SetExtent(cd->FrameCount(), box);

    return true;
}

bool OSPRayPKDGeometry::getDirtyCallback(core::Call& call) {
    auto os = dynamic_cast<CallOSPRayAPIObject*>(&call);
    auto cd = this->getDataSlot.CallAs<megamol::geocalls::MultiParticleDataCall>();

    if (cd == nullptr)
        return false;
    if (this->InterfaceIsDirtyNoReset()) {
        os->setDirty();
    }
    if (this->datahash != cd->DataHash()) {
        os->SetDataHash(cd->DataHash());
    }
    return true;
}


} // namespace megamol::ospray
