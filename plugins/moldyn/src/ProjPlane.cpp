#include "ProjPlane.h"

#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"
#include "mmcore/param/Vector3fParam.h"

#include "geometry_calls/MultiParticleDataCall.h"


megamol::moldyn::ProjPlane::ProjPlane()
        : out_data_slot_("outData", "")
        , in_data_slot_("inData", "")
        , plane_pos_("plane_pos", "")
        , plane_normal_("plane_normal", "")
        , sampling_width_("sampling_width", "")
        , temp_smooth_("temp_smooth", "") {
    in_data_slot_.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&in_data_slot_);

    plane_pos_ << new core::param::Vector3fParam(vislib::math::Vector<float, 3>(0, 50, 0));
    MakeSlotAvailable(&plane_pos_);

    plane_normal_ << new core::param::Vector3fParam(vislib::math::Vector<float, 3>(0, 1, 0));
    MakeSlotAvailable(&plane_normal_);

    sampling_width_ << new core::param::FloatParam(1.f, std::numeric_limits<float>::min());
    MakeSlotAvailable(&sampling_width_);

    temp_smooth_ << new core::param::IntParam(10, 1);
    MakeSlotAvailable(&temp_smooth_);
}


megamol::moldyn::ProjPlane::~ProjPlane() {
    this->Release();
}


bool megamol::moldyn::ProjPlane::create() {
    return true;
}


void megamol::moldyn::ProjPlane::release() {}


bool rayPlaneIntersection(
    glm::vec3 const& p_pos, glm::vec3 const& p_normal, glm::vec3 const& ray_org, glm::vec3 const& ray_dir, float& t) {
    auto const denom = glm::dot(p_normal, ray_dir);

    if (std::abs(denom) > 0.00001f) {
        t = glm::dot(p_pos - ray_org, p_normal) / denom;
        if (t >= 0.f)
            return true;
    }

    return false;
}


bool megamol::moldyn::ProjPlane::Render(core::view::CallRender3D& call) {
    auto in_data = in_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (in_data == nullptr)
        return false;

    if (!(*in_data)(0))
        return false;

    if (in_data->FrameID() != frame_id_ || in_data->DataHash() != in_data_hash_) {
        auto const pl_count = in_data->GetParticleListCount();

        size_t totalParts = 0;
        for (unsigned int i = 0; i < pl_count; ++i) {
            totalParts += in_data->AccessParticles(i).GetCount();
        }

        all_parts_.clear();
        all_parts_.reserve(totalParts);

        size_t allpartcnt = 0;
        for (unsigned int pli = 0; pli < pl_count; ++pli) {
            auto const& pl = in_data->AccessParticles(pli);
            /*if (!isListOK(in, pli) || !isDirOK(static_cast<metricsEnum>(theMetrics), in, pli)) {
                megamol::core::utility::log::Log::DefaultLog.WriteWarn(
                    "ParticleThermodyn: ignoring list %d because it either has no proper positions or no velocity",
                    pli);
                continue;
            }*/

            // unsigned int vert_stride = 0;
            // if (pl.GetVertexDataType() == DirectionalParticleDataCall::Particles::VertexDataType::VERTDATA_FLOAT_XYZ)
            // vert_stride = 12; else if (pl.GetVertexDataType() ==
            // DirectionalParticleDataCall::Particles::VertexDataType::VERTDATA_FLOAT_XYZR) vert_stride = 16; else
            // continue; vert_stride = std::max<unsigned int>(vert_stride, pl.GetVertexDataStride()); const unsigned char
            // *vert = static_cast<const unsigned char*>(pl.GetVertexData());

            auto const part_cnt = pl.GetCount();

            for (int part_i = 0; part_i < part_cnt; ++part_i) {
                all_parts_.push_back(allpartcnt + part_i);
            }
            allpartcnt += pl.GetCount();
        }

        myPts_ = std::make_shared<datatools::simplePointcloud>(in_data, all_parts_);

        particleTree_ = std::make_shared<my_kd_tree_t>(
            3 /* dim */, *myPts_, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
        particleTree_->buildIndex();

        // generate sampling positions:
        //   camera ray intersection with plane
        //   sample at these positions
        //   don't forget temporal averaging

        auto fbo = call.GetFramebuffer();

        auto cam = call.GetCamera();
        auto pose = cam.getPose();
        auto fov = cam.get<core::view::Camera::FieldOfViewY>();
        auto aspect = cam.get<core::view::Camera::AspectRatio>();
        auto tile = cam.get<core::view::Camera::ImagePlaneTile>();
        auto wh = glm::vec2(fbo->getWidth(), fbo->getHeight());
        //auto wh = tile.tile_end - tile.tile_start;
        auto num_pixels = wh.x * wh.y;
        std::vector<std::pair<glm::vec3, glm::vec3>> rays;
        rays.reserve(num_pixels);

        auto const near_clip = 100.f;
        auto const th = std::tan(0.5f * fov) * near_clip;
        auto const rw = th * aspect;


        for (unsigned int y = 0; y < wh.y; ++y) {
            for (unsigned int x = 0; x < wh.x; ++x) {
                float u = -rw + (rw + rw) * float(x) / wh.x;
                float v = -(th + (-th - th) * float(y) / wh.y);
                const glm::vec3 direction = (near_clip * pose.direction + u * pose.right + v * pose.up);
                rays.push_back(std::make_pair(pose.position, glm::normalize(direction)));
            }
        }


        auto vislib_pos = plane_pos_.Param<core::param::Vector3fParam>()->Value();
        auto vislib_normal = plane_normal_.Param<core::param::Vector3fParam>()->Value();

        auto const plane_pos = glm::vec3(vislib_pos.GetX(), vislib_pos.GetY(), vislib_pos.GetZ());
        auto const plane_normal = glm::vec3(vislib_normal.GetX(), vislib_normal.GetY(), vislib_normal.GetZ());

        std::vector<glm::vec3> sample_pos;
        sample_pos.reserve(num_pixels);

        for (auto const& [pos, dir] : rays) {
            float t = -1;
            rayPlaneIntersection(plane_pos, plane_normal, pos, dir, t);
            auto inter = pos + dir * t;
            sample_pos.push_back(inter);
        }

        auto const sample_width = sampling_width_.Param<core::param::FloatParam>()->Value();
        auto const sqr_sample_width = sample_width * sample_width;

        nanoflann::SearchParams params;
        params.sorted = false;

        std::vector<std::pair<size_t, float>> matches;
        matches.reserve(100);

        std::vector<glm::vec3> dir_samples;
        dir_samples.reserve(num_pixels);

        glm::vec3 min_vel(std::numeric_limits<float>::max());
        glm::vec3 max_vel(std::numeric_limits<float>::lowest());

        for (auto const& query : sample_pos) {
            /*std::vector<size_t> ind(10);
            std::vector<size_t> ind(10);*/
            auto N = particleTree_->radiusSearch(glm::value_ptr(query), sqr_sample_width, matches, params);
            //particleTree_->knnSearch(glm::value_ptr(query), 10, matches, params)
            glm::vec3 dir = glm::vec3(0);
            for (auto const& m : matches) {
                auto velo = myPts_->get_velocity(m.first);
                dir += glm::vec3(velo[0], velo[1], velo[2]);
            }
            if (!matches.empty())
                dir /= static_cast<float>(matches.size());
            dir_samples.push_back(dir);

            if (min_vel.x > dir.x)
                min_vel.x = dir.x;
            if (min_vel.y > dir.y)
                min_vel.y = dir.y;
            if (min_vel.z > dir.z)
                min_vel.z = dir.z;

            if (max_vel.x < dir.x)
                max_vel.x = dir.x;
            if (max_vel.y < dir.y)
                max_vel.y = dir.y;
            if (max_vel.z < dir.z)
                max_vel.z = dir.z;
        }

        auto vel_diff = max_vel - min_vel;

        std::vector<unsigned int> colors(num_pixels);
        for (unsigned int pix = 0; pix < num_pixels; ++pix) {
            auto const& dir = dir_samples[pix];
            auto const vel = (dir - min_vel) / vel_diff;
            colors[pix] = glm::packUint4x8(glm::u8vec4(vel * 255.f, 255));
        }

        fbo->colorBuffer = colors;

        std::vector<float> depths(num_pixels,0.5f);

        fbo->depthBuffer = depths;

        fbo->width = wh.x;
        fbo->height = wh.y;

        in_data_hash_ = in_data->DataHash();
        frame_id_ = in_data->FrameID();
    }

    return true;
}


bool megamol::moldyn::ProjPlane::GetExtents(core::view::CallRender3D& call) {
    auto in_data = in_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (in_data == nullptr)
        return false;

    if (!(*in_data)(1))
        return false;

    call.SetTimeFramesCount(in_data->FrameCount());

    call.AccessBoundingBoxes().SetBoundingBox(in_data->AccessBoundingBoxes().ObjectSpaceBBox());
    call.AccessBoundingBoxes().SetClipBox(in_data->AccessBoundingBoxes().ObjectSpaceClipBox());

    return true;
}


//bool megamol::moldyn::ProjPlane::get_data_cb(core::Call& c) {
//    auto in_data = in_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
//    if (in_data == nullptr)
//        return false;
//
//    if (!(*in_data)(0))
//        return false;
//
//    if (in_data->FrameID() != frame_id_ || in_data->DataHash() != in_data_hash_) {
//        auto const pl_count = in_data->GetParticleListCount();
//
//        size_t totalParts = 0;
//        for (unsigned int i = 0; i < pl_count; ++i) {
//            totalParts += in_data->AccessParticles(i).GetCount();
//        }
//
//        all_parts_.clear();
//        all_parts_.reserve(totalParts);
//
//        size_t allpartcnt = 0;
//        for (unsigned int pli = 0; pli < pl_count; ++pli) {
//            auto const& pl = in_data->AccessParticles(pli);
//            /*if (!isListOK(in, pli) || !isDirOK(static_cast<metricsEnum>(theMetrics), in, pli)) {
//                megamol::core::utility::log::Log::DefaultLog.WriteWarn(
//                    "ParticleThermodyn: ignoring list %d because it either has no proper positions or no velocity",
//                    pli);
//                continue;
//            }*/
//
//            // unsigned int vert_stride = 0;
//            // if (pl.GetVertexDataType() == DirectionalParticleDataCall::Particles::VertexDataType::VERTDATA_FLOAT_XYZ)
//            // vert_stride = 12; else if (pl.GetVertexDataType() ==
//            // DirectionalParticleDataCall::Particles::VertexDataType::VERTDATA_FLOAT_XYZR) vert_stride = 16; else
//            // continue; vert_stride = std::max<unsigned int>(vert_stride, pl.GetVertexDataStride()); const unsigned char
//            // *vert = static_cast<const unsigned char*>(pl.GetVertexData());
//
//            auto const part_cnt = pl.GetCount();
//
//            for (int part_i = 0; part_i < part_cnt; ++part_i) {
//                all_parts_.push_back(allpartcnt + part_i);
//            }
//            allpartcnt += pl.GetCount();
//        }
//
//        myPts_ = std::make_shared<datatools::simplePointcloud>(*in_data, all_parts_);
//
//        particleTree_ = std::make_shared<my_kd_tree_t>(
//            3 /* dim */, *myPts_, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
//        particleTree_->buildIndex();
//
//        // generate sampling positions:
//        //   camera ray intersection with plane
//        //   sample at these positions
//        //   don't forget temporal averaging
//
//        in_data_hash_ = in_data->DataHash();
//        frame_id_ = in_data->FrameID();
//    }
//
//    return true;
//}
//
//
//bool megamol::moldyn::ProjPlane::get_extent_cb(core::Call& c) {
//    return true;
//}
