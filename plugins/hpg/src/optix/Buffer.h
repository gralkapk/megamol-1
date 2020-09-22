#pragma once

#include <algorithm>
#include <cstddef>

#include "cuda.h"
#include "cudaGL.h"

namespace megamol::hpg::optix {

class DeviceBuffer {
public:
    DeviceBuffer() : _size(0), _ptr(0) {}

    DeviceBuffer(std::size_t size) : _size(size) { cuMemAlloc(&get_pointer(), size); }

    void write(void* src, std::size_t size) { cuMemcpyHtoD(_ptr, src, size); }

    void write_async(void* src, std::size_t size, CUstream stream) { cuMemcpyHtoDAsync(_ptr, src, size, stream); }

    void read(void* dst, std::size_t size) const { cuMemcpyDtoH(dst, _ptr, size); }

    void read_async(void* dst, std::size_t size, CUstream stream) const { cuMemcpyDtoHAsync(dst, _ptr, size, stream); }

    virtual ~DeviceBuffer() { cuMemFree(get_pointer()); };

protected:
    std::size_t& get_size() { return _size; }

    CUdeviceptr& get_pointer() { return _ptr; }

private:
    std::size_t _size;

    CUdeviceptr _ptr;
}; // end class DeviceBuffer


class PinnedBuffer : public DeviceBuffer {
public:
    PinnedBuffer(std::size_t size) {
        cuMemHostAlloc(&_host_ptr, size, CU_MEMHOSTALLOC_WRITECOMBINED | CU_MEMHOSTALLOC_DEVICEMAP);
        cuMemHostGetDevicePointer(&get_pointer(), _host_ptr, 0);
        get_size() = size;
    }

    virtual ~PinnedBuffer() {
        cuMemFreeHost(_host_ptr);
        get_pointer() = 0;
    }

private:
    void* _host_ptr;
}; // end class PinnedBuffer


class GraphicsBuffer : public DeviceBuffer {
public:
    void map(CUstream stream) { cuGraphicsMapResources(1, &_res, stream); }

    void unmap(CUstream stream) { cuGraphicsUnmapResources(1, &_res, stream); }

    virtual ~GraphicsBuffer() { cuGraphicsUnregisterResource(_res); }

protected:
    CUgraphicsResource& get_resource() { return _res; }

    void get_mapped_pointer() { cuGraphicsResourceGetMappedPointer(&get_pointer(), &get_size(), _res); }

private:
    CUgraphicsResource _res;
}; // end class GraphicsBuffer


class GraphicsTextureBuffer : public GraphicsBuffer {
public:
    GraphicsTextureBuffer(GLuint image, GLenum target, unsigned int flags) {
        cuGraphicsGLRegisterImage(&get_resource(), image, target, flags);
        get_mapped_pointer();
    }

private:
}; // end class GraphicsTextureBuffer


class GraphicsBOBuffer : public GraphicsBuffer {
public:
    GraphicsBOBuffer(GLuint buffer, unsigned int flags) {
        cuGraphicsGLRegisterBuffer(&get_resource(), buffer, flags);
        get_mapped_pointer();
    }

private:
}; // end class GraphicsBOBuffer


} // end namespace megamol::hpg::optix