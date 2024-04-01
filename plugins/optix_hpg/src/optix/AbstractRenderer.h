#pragma once

#include "mmstd/renderer/RendererModule.h"

#include "CallRender3DCUDA.h"

#include "optix/Utils.h"

#include "framestate.h"

#include "CUDA_Context.h"
#include "SBT.h"

#include "optix/Context.h"

#include "optix/CallGeometry.h"

namespace megamol::optix_hpg {
template<typename FSType>
class AbstractRenderer : public core::view::RendererModule<CallRender3DCUDA, core::Module> {
public:
    static void requested_lifetime_resources(frontend_resources::ResourceRequest& req) {
        core::view::RendererModule<CallRender3DCUDA, Module>::requested_lifetime_resources(req);
        req.require<frontend_resources::CUDA_Context>();
    }

    AbstractRenderer();

    virtual ~AbstractRenderer();

    bool Render(CallRender3DCUDA& call) override;

    bool GetExtents(CallRender3DCUDA& call) override;

protected:
    bool create() override;

    void release() override;

    FSType& get_framestate() {
        return frame_state_;
    }

    std::unique_ptr<Context> const& get_context() const {
        return optix_ctx_;
    }

    MMOptixSBT& get_sbt() {
        return sbt_;
    }

    OptixPipeline& get_pipeline() {
        return pipeline_;
    }

    CUdeviceptr const& get_frame_state_buffer() const {
        return frame_state_buffer_;
    }

private:
    virtual void setup() = 0;

    virtual bool is_dirty() = 0;

    virtual void reset_dirty() = 0;

    virtual void on_cam_param_change(
        core::view::Camera const& cam, core::view::Camera::PerspectiveParameters const& cam_intrinsics) = 0;

    virtual void on_cam_pose_change(core::view::Camera::Pose const& cam_pose) = 0;

    virtual void on_change_data(OptixTraversableHandle world) = 0;

    virtual void on_change_background(glm::vec4 const& bg) = 0;

    virtual void on_change_programs(std::tuple<OptixProgramGroup const*, uint32_t> const& programs) = 0;

    virtual void on_change_parameters() = 0;

    virtual void on_change_viewport(glm::uvec2 const& viewport, std::shared_ptr<CUDAFramebuffer> fbo) = 0;

    virtual void on_change_sbt(std::tuple<void const*, uint32_t, uint64_t> const& records) = 0;

    virtual void on_launch_func(){};

    core::CallerSlot in_geo_slot_;

    std::unique_ptr<Context> optix_ctx_;

    CUdeviceptr frame_state_buffer_;

    MMOptixSBT sbt_;

    OptixPipeline pipeline_ = 0;

    FSType frame_state_;

    glm::uvec2 current_fb_size_;

    unsigned int frame_id_ = std::numeric_limits<unsigned int>::max();

    std::size_t in_data_hash_ = std::numeric_limits<std::size_t>::max();

    core::view::Camera::Pose old_cam_pose_;

    core::view::Camera::PerspectiveParameters old_cam_intrinsics_;

    glm::vec4 old_bg_ = glm::vec4(-1);
};

template<typename FSType>
megamol::optix_hpg::AbstractRenderer<FSType>::AbstractRenderer() : in_geo_slot_("inGeo", "") {
    in_geo_slot_.SetCompatibleCall<CallGeometryDescription>();
    MakeSlotAvailable(&in_geo_slot_);

    // Callback should already be set by RendererModule
    this->MakeSlotAvailable(&this->chainRenderSlot);

    // Callback should already be set by RendererModule
    this->MakeSlotAvailable(&this->renderSlot);
}

template<typename FSType>
megamol::optix_hpg::AbstractRenderer<FSType>::~AbstractRenderer() {
    this->Release();
}

template<typename FSType>
bool megamol::optix_hpg::AbstractRenderer<FSType>::Render(CallRender3DCUDA& call) {
    auto const viewport = glm::uvec2(call.GetFramebuffer()->width, call.GetFramebuffer()->height);

    call.GetFramebuffer()->data.exec_stream = optix_ctx_->GetExecStream();

    auto in_geo = in_geo_slot_.CallAs<CallGeometry>();
    if (in_geo == nullptr)
        return false;

    in_geo->set_ctx(optix_ctx_.get());
    if (!(*in_geo)())
        return false;

    bool rebuild_sbt = false;

    // change viewport
    if (viewport != current_fb_size_) {
        on_change_viewport(viewport, call.GetFramebuffer());

        rebuild_sbt = true;

        current_fb_size_ = viewport;
    }

    // Camera
    core::view::Camera cam = call.GetCamera();

    auto const cam_pose = cam.get<core::view::Camera::Pose>();
    auto const cam_intrinsics = cam.get<core::view::Camera::PerspectiveParameters>();

    // change camera parameters
    if (!(cam_intrinsics == old_cam_intrinsics_)) {
        on_cam_param_change(cam, cam_intrinsics);
        old_cam_intrinsics_ = cam_intrinsics;
    }

    // change camera pose
    if (!(cam_pose == old_cam_pose_)) {
        on_cam_pose_change(cam_pose);
        old_cam_pose_ = cam_pose;
    }

    // update parameters
    if (is_dirty()) {
        on_change_parameters();
        reset_dirty();
    }

    // change background color
    if (old_bg_ != call.BackgroundColor()) {
        on_change_background(call.BackgroundColor());
        rebuild_sbt = true;
        old_bg_ = call.BackgroundColor();
    }

    // change data
    if (frame_id_ != in_geo->FrameID() || in_data_hash_ != in_geo->DataHash() || in_geo->has_geo_update()) {
        on_change_data(*in_geo->get_handle());

        rebuild_sbt = true;
        frame_id_ = in_geo->FrameID();
        in_data_hash_ = in_geo->DataHash();
    }

    // change programs
    if (in_geo->has_program_update()) {
        on_change_programs(in_geo->get_program_groups());
    }

    CUDA_CHECK_ERROR(
        cuMemcpyHtoDAsync(frame_state_buffer_, &frame_state_, sizeof(frame_state_), optix_ctx_->GetExecStream()));

    if (rebuild_sbt || in_geo->has_sbt_update()) {
        on_change_sbt(in_geo->get_record());
    }

    OPTIX_CHECK_ERROR(optixLaunch(pipeline_, optix_ctx_->GetExecStream(), 0, 0, sbt_, viewport.x, viewport.y, 1));

    // post launch callback
    on_launch_func();

    ++frame_state_.frameIdx;

    return true;
}

template<typename FSType>
bool megamol::optix_hpg::AbstractRenderer<FSType>::GetExtents(CallRender3DCUDA& call) {
    auto in_geo = in_geo_slot_.CallAs<CallGeometry>();
    if (in_geo != nullptr) {
        in_geo->SetFrameID(static_cast<unsigned int>(call.Time()));
        if (!(*in_geo)(1))
            return false;
        call.SetTimeFramesCount(in_geo->FrameCount());

        call.AccessBoundingBoxes() = in_geo->AccessBoundingBoxes();
    } else {
        call.SetTimeFramesCount(1);
        call.AccessBoundingBoxes().Clear();
    }

    return true;
}

template<typename FSType>
bool megamol::optix_hpg::AbstractRenderer<FSType>::create() {
    auto& cuda_res = frontend_resources.get<frontend_resources::CUDA_Context>();
    if (cuda_res.ctx_ != nullptr) {
        optix_ctx_ = std::make_unique<Context>(cuda_res);
    } else {
        return false;
    }

    on_change_parameters();

    CUDA_CHECK_ERROR(cuMemAlloc(&frame_state_buffer_, sizeof(device::FrameState)));

    setup();

    return true;
}

template<typename FSType>
void megamol::optix_hpg::AbstractRenderer<FSType>::release() {
    CUDA_CHECK_ERROR(cuMemFree(frame_state_buffer_));
    OPTIX_CHECK_ERROR(optixPipelineDestroy(pipeline_));
}
} // namespace megamol::optix_hpg
