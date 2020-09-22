#pragma once

#include "optix.h"
#include "optix_stubs.h"

#include "mmcore/utility/log/Log.h"

namespace megamol::hpg::optix {

inline void optix_log_cb(unsigned int level, const char* tag, const char* message, void* cbdata) {
    megamol::core::utility::log::Log::DefaultLog.WriteMsg(level, "OptiX Msg with tag: %s; Msg: %s", tag, message);
}

inline CUresult cu_print_error(CUresult error) {
    if (error != CUDA_SUCCESS) {
        char const* error_name = nullptr;
        cuGetErrorName(error, &error_name);
        char const* error_string = nullptr;
        cuGetErrorString(error, &error_string);
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "%s: %s in %s:%d", error_name, error_string, __FILE__, __LINE__);
    }
    return error;
}

class OptiXContext {
public:
    OptiXContext();
    ~OptiXContext();

    operator OptixDeviceContext() { return _optix_ctx; }

private:
    CUstream _stream;

    OptixDeviceContext _optix_ctx = 0;
};

} // end namespace megamol::hpg::optix