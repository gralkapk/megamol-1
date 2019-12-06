#pragma once

#include <array>
#include <stdexcept>
#include <vector>
#include "cuda.h"
#include "cuda_runtime.h"

namespace megamol {
namespace hpg {
namespace optix {

template <typename T> class Buffer {
public:
    Buffer(size_t size) : ptr_(nullptr), size_(size) { cudaMalloc(&ptr_, size_); }

    Buffer(Buffer const& rhs) = delete;

    Buffer& operator=(Buffer const& rhs) = delete;

    Buffer(Buffer&& rhs) noexcept : ptr_(nullptr), size_(0) { *this = std::move(rhs); }

    Buffer& operator=(Buffer&& rhs) noexcept {
        if (std::addressof(rhs) != this) {
            std::swap(ptr_, rhs.ptr_);
            std::swap(size_, rhs.size_);
        }
        return *this;
    }

    ~Buffer() { cudaFree(ptr_); }

    bool Upload(T const* ptr, size_t size) {
        if (size <= size_) {
            return cudaMemcpy(ptr_, ptr, size, cudaMemcpyHostToDevice) == cudaSuccess;
        }
        return false;
    }

    bool Download(T* ptr, size_t size) {
        if (size <= size_) {
            return cudaMemcpy(ptr, ptr_, size, cudaMemcpyDeviceToHost);
        }
        return false;
    }

    bool Upload(std::vector<T> const& data) { return Upload(data.data(), data.size()); }

    template <size_t N> bool Upload(std::array<T, N> const& data) { return Upload(data.data(), N); }

    bool Download(std::vector<T>& data) { return Download(data.data(), data.size()); }

    template <size_t N> bool Download(std::array<T, N>& data) { return Download(data.data(), N); }

    [[nodiscard]] size_t Size() const { return size_; }

    operator CUdeviceptr() const { return reinterpret_cast<CUdeviceptr>(ptr_); }

    operator void*() const { return reinterpret_cast<void*>(ptr_); }

private:
    T* ptr_;
    size_t size_;
};

} // end namespace optix
} // end namespace hpg
} // end namespace megamol