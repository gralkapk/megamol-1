#pragma once

#include <memory>
#include "cuda.h"
#include "cuda_runtime.h"
#include "optix.h"

namespace megamol {
namespace hpg {

class OptixContext {
public:
    OptixContext() : context_(nullptr) {
        cudaFree(nullptr);
        CUcontext cuCtx = nullptr;
        optixDeviceContextCreate(cuCtx, nullptr, &context_);
    }

    OptixContext(OptixContext const& rhs) = delete;

    OptixContext& operator=(OptixContext& rhs) = delete;

    OptixContext(OptixContext&& rhs) noexcept : context_(nullptr) { *this = std::move(rhs); }

    OptixContext& operator=(OptixContext&& rhs) noexcept {
        if (std::addressof(rhs) != this) {
            std::swap(context_, rhs.context_);
        }
        return *this;
    }

    ~OptixContext() { optixDeviceContextDestroy(context_); }

    operator OptixDeviceContext() const { return this->context_; }

private:
    OptixDeviceContext context_;
};


template <typename C> class Context {
public:
private:
    C context_;
};


} // end namespace hpg
} // end namespace megamol