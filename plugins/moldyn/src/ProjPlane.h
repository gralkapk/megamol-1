#pragma once

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "mmstd/renderer/Renderer3DModule.h"

#include "datatools/PointcloudHelpers.h"
#include "nanoflann.hpp"

namespace megamol::moldyn {
class ProjPlane : public core::view::Renderer3DModule {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() {
        return "ProjPlane";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() {
        return "ProjPlane";
    }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() {
        return true;
    }

    ProjPlane();

    virtual ~ProjPlane();

protected:
    bool create() override;

    void release() override;

private:
    bool Render(core::view::CallRender3D& call) override;

    bool GetExtents(core::view::CallRender3D& call) override;

    /*bool get_data_cb(core::Call& c);

    bool get_extent_cb(core::Call& c);*/

    core::CalleeSlot out_data_slot_;

    core::CallerSlot in_data_slot_;

    core::param::ParamSlot plane_pos_;

    core::param::ParamSlot plane_normal_;

    core::param::ParamSlot sampling_width_;

    core::param::ParamSlot temp_smooth_;

    unsigned int frame_id_ = 0;

    uint64_t in_data_hash_ = std::numeric_limits<uint64_t>::max();

    uint64_t out_data_hash_ = 0;

    typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, datatools::simplePointcloud>,
        datatools::simplePointcloud, 3 /* dim */, std::size_t>
        my_kd_tree_t;

    std::shared_ptr<my_kd_tree_t> particleTree_;
    std::shared_ptr<datatools::simplePointcloud> myPts_;

    std::vector<size_t> all_parts_;
};
} // namespace megamol::moldyn
