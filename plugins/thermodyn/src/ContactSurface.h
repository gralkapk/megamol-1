#pragma once

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"

#include <nanoflann.hpp>

#include "datatools/PointcloudHelpers.h"

namespace megamol::thermodyn {
class ContactSurface : public core::Module {
public:
    /** Return module class name */
    static const char* ClassName(void) {
        return "ContactSurface";
    }

    /** Return module class description */
    static const char* Description(void) {
        return "ContactSurface";
    }

    /** Module is always available */
    static bool IsAvailable(void) {
        return true;
    }

    ContactSurface();

    virtual ~ContactSurface();

protected:
    bool create() override;

    void release() override;

private:
    core::CalleeSlot out_data_slot_;

    core::CallerSlot in_data_slot_;

    core::param::ParamSlot knn_slot_;

    core::param::ParamSlot distance_threshold_slot_;

    bool get_data_cb(core::Call& c);

    bool get_extent_cb(core::Call& c);

    bool check_dirty() const {
        return knn_slot_.IsDirty() || distance_threshold_slot_.IsDirty();
    }

    void reset_dirty() {
        knn_slot_.ResetDirty();
        distance_threshold_slot_.ResetDirty();
    }

    uint64_t in_data_hash_ = std::numeric_limits<uint64_t>::max();

    uint64_t out_data_hash_ = 0;

    unsigned int frame_id_ = -1;

    typedef datatools::genericPointcloud<float, 3> pointcloud_t;

    typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, pointcloud_t>, pointcloud_t,
        3 /* dim */, std::size_t>
        my_kd_tree_t;

    std::shared_ptr<my_kd_tree_t> particleTree;
    std::shared_ptr<pointcloud_t> myPts;

    std::vector<float> out_data_vec_;
};
} // namespace megamol::thermodyn
