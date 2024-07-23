#include "CUDAQuery.h"

#ifdef MEGAMOL_USE_CUDA
#include <cuda.h>
#endif

namespace megamol::frontend_resources::performance {
CUDAQuery::CUDAQuery() {
#ifdef MEGAMOL_USE_CUDA
    
#endif
}

CUDAQuery::~CUDAQuery() {
#ifdef MEGAMOL_USE_CUDA
    
#endif
}

void CUDAQuery::Counter(void* userData) {
#ifdef MEGAMOL_USE_CUDA
    cuLaunchHostFunc(
        reinterpret_cast<CUstream>(userData),
        [](void* v) {
            *reinterpret_cast<std::chrono::steady_clock::time_point*>(v) = std::chrono::steady_clock::now();
        },
        &value_);
#endif
}

time_point CUDAQuery::GetNW() {
#ifdef MEGAMOL_USE_OPENGL
    
#endif
    return value_;
}
} // namespace megamol::frontend_resources::performance
