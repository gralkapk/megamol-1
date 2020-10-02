#include "stdafx.h"
#include "Pipeline.h"

#include "hpg/optix/CallContext.h"
#include "hpg/optix/CallProgramGroup.h"
#include "hpg/optix/CallPipeline.h"

#include "OptiXCompileOptions.h"


megamol::hpg::optix::Pipeline::Pipeline()
    : _out_pipeline_slot("outPipe", ""), _in_context_slot("inContext", ""), _in_groups_slot("inGroups", "") {
    _out_pipeline_slot.SetCallback(
        CallPipeline::ClassName(), CallPipeline::FunctionName(0), &Pipeline::get_pipeline_cb);
    MakeSlotAvailable(&_out_pipeline_slot);

    _in_context_slot.SetCompatibleCall<CallContextDescription>();
    MakeSlotAvailable(&_in_context_slot);

    _in_groups_slot.SetCompatibleCall<CallProgramGroupDescription>();
    MakeSlotAvailable(&_in_groups_slot);

    _pipe = std::make_shared<OptixPipeline>(nullptr);
}


megamol::hpg::optix::Pipeline::~Pipeline() { this->Release(); }


bool megamol::hpg::optix::Pipeline::create() { return true; }


void megamol::hpg::optix::Pipeline::release() {}


bool megamol::hpg::optix::Pipeline::get_pipeline_cb(core::Call& c) {
    auto out_call = dynamic_cast<CallPipeline*>(&c);
    if (out_call == nullptr) return false;
    auto in_context = _in_context_slot.CallAs<CallContext>();
    if (in_context == nullptr) return false;
    auto in_groups = _in_groups_slot.CallAs<CallProgramGroup>();
    if (in_groups == nullptr) return false;

    if (in_context->is_dirty() || in_groups->is_dirty()) {
        auto ctx = in_context->get_context();
        if (ctx == nullptr) return false;
        auto groups = in_groups->get_groups();
        if (groups == nullptr) return false;

        auto options = OptiXCompileOptions();

        optixPipelineCreate(*ctx, &static_cast<OptixPipelineCompileOptions>(options),
            &static_cast<OptixPipelineLinkOptions>(options), groups->data(), groups->size(), nullptr, nullptr,
            _pipe.get());

        in_groups->reset_dirty();
        in_context->reset_dirty();

        out_call->set_dirty();
    }

    out_call->set_pipeline(_pipe);

    return true;
}
