#pragma once

#include "mmcore/utility/log/Log.h"

#include "cuda.h"

#include "optix.h"
#include "optix_stubs.h"

namespace megamol::hpg::optix {
inline CUresult print_cuda_error(CUresult ec, char const* file, int line, const char* func) {
    if (ec != CUDA_SUCCESS) {
        const char* en = nullptr;
        cuGetErrorName(ec, &en);
        const char* es = nullptr;
        cuGetErrorString(ec, &es);
        megamol::core::utility::log::Log::DefaultLog.WriteError("[%s] at %s:%d ... (%s) %s", func, file, line, en, es);
    }
    return ec;
}

inline OptixResult print_optix_error(OptixResult ec, char const* file, int line, const char* func) {
    if (ec != OPTIX_SUCCESS) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[%s] at %s:%d ... (%s) %s", func, file, line, optixGetErrorName(ec), optixGetErrorString(ec));
    }
    return ec;
}

inline void optix_log_callback(unsigned int level, char const* tag, char const* message, void* cbdata) {
    megamol::core::utility::log::Log::DefaultLog.WriteInfo("[OptiX Debug Msg %s] (%d) %s", tag, level, message);
}
} // namespace megamol::hpg::optix

#ifdef DEBUG
#define CUDA_CHECK_ERROR(x) print_cuda_error((x), __FILE__, __LINE__, __func__)
#define OPTIX_CHECK_ERROR(x) print_optix_error((x), __FILE__, __LINE__, __func__)
#else
#define CUDA_CHECK_ERROR(x) (x)
#define OPTIX_CHECK_ERROR(x) (x)
#endif
