#include "CUDAQuery.h"

#ifdef MEGAMOL_USE_CUDA
#include <cuda.h>
#endif

namespace megamol::frontend_resources::performance {
CUDAQuery::CUDAQuery() {
#ifdef MEGAMOL_USE_CUDA
    cuEventCreate(reinterpret_cast<CUevent*>(&handle_), CU_EVENT_DEFAULT);
#endif
}

CUDAQuery::~CUDAQuery() {
#ifdef MEGAMOL_USE_CUDA
    cuEventDestroy(reinterpret_cast<CUevent>(handle_));
#endif
}

void CUDAQuery::Counter(void* userData) {
#ifdef MEGAMOL_USE_CUDA
    cuEventRecord(reinterpret_cast<CUevent>(handle_), reinterpret_cast<CUstream>(userData));
#endif
}

time_point CUDAQuery::GetNW() {
    return value_;
}

void CUDAQuery::Sync(std::shared_ptr<AnyQuery> start, void* userData) {
    cuStreamSynchronize(reinterpret_cast<CUstream>(userData));
    float time_in_ms = 0;
    cuEventElapsedTime(&time_in_ms, reinterpret_cast<CUevent>(start->GetHandle()), reinterpret_cast<CUevent>(handle_));
    start->SetValue(time_point{std::chrono::nanoseconds(0)});
    value_ = time_point{std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::duration<float>(time_in_ms/1000.f))};
}
} // namespace megamol::frontend_resources::performance
