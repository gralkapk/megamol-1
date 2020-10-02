#pragma once

#include "optix.h"

namespace megamol::hpg::optix {

class OptiXCompileOptions {
public:
    OptiXCompileOptions() : _module_options{0}, _pipeline_options{0}, _pipeline_link_options{0} {
        _module_options.optLevel = OptixCompileOptimizationLevel::OPTIX_COMPILE_OPTIMIZATION_LEVEL_3;
        _module_options.debugLevel = OptixCompileDebugLevel::OPTIX_COMPILE_DEBUG_LEVEL_NONE;

        _pipeline_options.usesMotionBlur = false;
        _pipeline_options.traversableGraphFlags = OptixTraversableGraphFlags::OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_ANY;
        _pipeline_options.numPayloadValues = 1;
        _pipeline_options.numAttributeValues = 2;
        _pipeline_options.exceptionFlags = OptixExceptionFlags::OPTIX_EXCEPTION_FLAG_NONE;
        _pipeline_options.pipelineLaunchParamsVariableName = "params";
        _pipeline_options.usesPrimitiveTypeFlags = OptixPrimitiveTypeFlags::OPTIX_PRIMITIVE_TYPE_FLAGS_CUSTOM;

        _pipeline_link_options.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_NONE;
        _pipeline_link_options.maxTraceDepth = 1;
    }

    void set_maxRegisterCount(int reg_count) { _module_options.maxRegisterCount = reg_count; }

    void set_optLevel(OptixCompileOptimizationLevel opt_lvl) { _module_options.optLevel = opt_lvl; }

    void set_debugLevel(OptixCompileDebugLevel dbg_lvl) {
        _module_options.debugLevel = dbg_lvl;
        _pipeline_link_options.debugLevel = dbg_lvl;
    }

    void set_usesMotionBlur(bool mo_blur) { _pipeline_options.usesMotionBlur = mo_blur; }

    void set_traversableGraphFlags(OptixTraversableGraphFlags tg_flags) {
        _pipeline_options.traversableGraphFlags = tg_flags;
    }

    void set_numPayloadValues(int num_plv) { _pipeline_options.numPayloadValues = num_plv; }

    void set_numAttributeValues(int num_av) { _pipeline_options.numAttributeValues = num_av; }

    void set_exceptionFlags(OptixExceptionFlags exp_flags) { _pipeline_options.exceptionFlags = exp_flags; }

    void set_pipelineLaunchParamsVariableName(char const* pv_name) {
        _pipeline_options.pipelineLaunchParamsVariableName = pv_name;
    }

    void set_usesPrimitiveTypeFlags(OptixPrimitiveTypeFlags pt_flags) {
        _pipeline_options.usesPrimitiveTypeFlags = pt_flags;
    }

    void set_maxTraceDepth(unsigned int depth) { _pipeline_link_options.maxTraceDepth = depth; }

    operator OptixModuleCompileOptions const &() const { return _module_options; }

    operator OptixPipelineCompileOptions const &() const { return _pipeline_options; }

    operator OptixPipelineLinkOptions const &() const { return _pipeline_link_options; }

private:
    OptixModuleCompileOptions _module_options;

    OptixPipelineCompileOptions _pipeline_options;

    OptixPipelineLinkOptions _pipeline_link_options;
};

} // end namespace megamol::hpg::optix