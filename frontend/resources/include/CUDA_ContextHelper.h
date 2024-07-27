#pragma once

#ifdef MEGAMOL_USE_CUDA
#include <cuda.h>
#include "mmcore/utility/log/Log.h"

inline CUresult print_cuda_error(CUresult ec, char const* file, int line) {
    if (ec != CUDA_SUCCESS) {
        const char* en = nullptr;
        cuGetErrorName(ec, &en);
        const char* es = nullptr;
        cuGetErrorString(ec, &es);
        megamol::core::utility::log::Log::DefaultLog.WriteError("CUDA Error at %s:%d ... (%s) %s", file, line, en, es);
    }
    return ec;
}

#ifdef DEBUG
#define CUDA_CHECK_ERROR(x) print_cuda_error((x), __FILE__, __LINE__)
#else
#define CUDA_CHECK_ERROR(x) (x)
#endif
#endif
