#include "AbstractRenderer.h"

#include "mmcore/param/BoolParam.h"

#include "optix/CallGeometry.h"

#ifdef MEGAMOL_USE_TRACY
#include <tracy/Tracy.hpp>
#endif


megamol::optix_hpg::AbstractRenderer::AbstractRenderer() : in_geo_slot_("inGeo", "")/*, profiling_slot_("profiling::enable", "")*/ {
    in_geo_slot_.SetCompatibleCall<CallGeometryDescription>();
    MakeSlotAvailable(&in_geo_slot_);

    // Callback should already be set by RendererModule
    this->MakeSlotAvailable(&this->chainRenderSlot);

    // Callback should already be set by RendererModule
    this->MakeSlotAvailable(&this->renderSlot);

    /*profiling_slot_ << new core::param::BoolParam(false);
    MakeSlotAvailable(&profiling_slot_);*/
}


megamol::optix_hpg::AbstractRenderer::~AbstractRenderer() {
    this->Release();
}


bool megamol::optix_hpg::AbstractRenderer::Render(CallRender3DCUDA& call) {
#ifdef MEGAMOL_USE_TRACY
    ZoneScopedN("AbstractRenderer::Render");
#endif
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

 #if 0
    if (profiling_slot_.Param<core::param::BoolParam>()->Value()) {
#ifdef MEGAMOL_USE_PROFILING
        auto& region = perf_man_->start_timer(launch_timer_, optix_ctx_->GetExecStream());
#endif
        CUDA_CHECK_ERROR(cuEventRecord(rend_start, optix_ctx_->GetExecStream()));
        OPTIX_CHECK_ERROR(optixLaunch(pipeline_, optix_ctx_->GetExecStream(), 0, 0, sbt_, viewport.x, viewport.y, 1));
        CUDA_CHECK_ERROR(cuEventRecord(rend_stop, optix_ctx_->GetExecStream()));
#ifdef MEGAMOL_USE_PROFILING
        region.end_region(optix_ctx_->GetExecStream());
#endif

        CUDA_CHECK_ERROR(cuStreamSynchronize(optix_ctx_->GetExecStream()));
        //CUDA_CHECK_ERROR(cuStreamWaitEvent(optix_ctx_->GetExecStream(), rend_start, CU_EVENT_WAIT_DEFAULT));
        //CUDA_CHECK_ERROR(cuStreamWaitEvent(optix_ctx_->GetExecStream(), rend_stop, CU_EVENT_WAIT_DEFAULT));
        float time_in_ms = 0.f;
        CUDA_CHECK_ERROR(cuEventElapsedTime(&time_in_ms, rend_start, rend_stop));
        core::utility::log::Log::DefaultLog.WriteInfo("[OptiX] Render time: %f ms", time_in_ms);
    } else {
        OPTIX_CHECK_ERROR(optixLaunch(pipeline_, optix_ctx_->GetExecStream(), 0, 0, sbt_, viewport.x, viewport.y, 1));
    }
#endif

#ifdef MEGAMOL_USE_PROFILING
    auto& region = perf_man_->start_timer(launch_timer_, optix_ctx_->GetExecStream());
#endif
    OPTIX_CHECK_ERROR(optixLaunch(pipeline_, optix_ctx_->GetExecStream(), 0, 0, sbt_, viewport.x, viewport.y, 1));
#ifdef MEGAMOL_USE_PROFILING
    region.end_region(optix_ctx_->GetExecStream());
#endif

    ++frame_state_.frameIdx;

    return true;
}


bool megamol::optix_hpg::AbstractRenderer::GetExtents(CallRender3DCUDA& call) {
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


bool megamol::optix_hpg::AbstractRenderer::create() {
    auto& cuda_res = frontend_resources.get<frontend_resources::CUDA_Context>();
    if (cuda_res.ctx_ != nullptr) {
        optix_ctx_ = std::make_unique<Context>(cuda_res);
    } else {
        return false;
    }

    on_change_parameters();

    CUDA_CHECK_ERROR(cuMemAlloc(&frame_state_buffer_, sizeof(device::FrameState)));

    setup();

    /*CUDA_CHECK_ERROR(cuEventCreate(&rend_start, CU_EVENT_BLOCKING_SYNC));
    CUDA_CHECK_ERROR(cuEventCreate(&rend_stop, CU_EVENT_BLOCKING_SYNC));*/

#ifdef MEGAMOL_USE_PROFILING
    perf_man_ = &const_cast<frontend_resources::performance::PerformanceManager&>(
        frontend_resources.get<frontend_resources::performance::PerformanceManager>());
    launch_timer_ = perf_man_->add_timer("optixLaunch", frontend_resources::performance::query_api::CUDA);
#endif

    return true;
}


void megamol::optix_hpg::AbstractRenderer::release() {
    CUDA_CHECK_ERROR(cuMemFree(frame_state_buffer_));
    OPTIX_CHECK_ERROR(optixPipelineDestroy(pipeline_));

    /*CUDA_CHECK_ERROR(cuEventDestroy(rend_start));
    CUDA_CHECK_ERROR(cuEventDestroy(rend_stop));*/
}
