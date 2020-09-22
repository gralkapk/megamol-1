#include "stdafx.h"
#include "optix/OptiXContext.h"


megamol::hpg::optix::OptiXContext::OptiXContext() {
    cu_print_error(cuInit(0));

    cu_print_error(cuStreamCreate(&_stream, CU_STREAM_NON_BLOCKING));

    OptixDeviceContextOptions opt = {optix_log_cb, nullptr, 3};

    optixDeviceContextCreate(0, &opt, &_optix_ctx);
}


megamol::hpg::optix::OptiXContext::~OptiXContext() {
    cuStreamDestroy(_stream);
    optixDeviceContextDestroy(_optix_ctx);
}
