#pragma once

#ifndef __CUDACC__
#define CU_CALLABLE
#else
#define CU_CALLABLE __host__ __device__
#endif

namespace megamol {
namespace optix_hpg {
namespace device {
struct BTParticle {
    unsigned int dim : 2;
    unsigned int x : 10;
    unsigned int y : 10;
    unsigned int z : 10;

    CU_CALLABLE unsigned int get(unsigned int dim) {
        if (dim == 0) {
            return x;
        }
        if (dim == 1) {
            return y;
        }
        if (dim == 2) {
            return z;
        }
    }

    CU_CALLABLE void to(glm::vec3 const& pos, glm::vec3 const& span, glm::vec3 const& lower) {
        auto const diff = (pos - lower) / span;
        x = static_cast<unsigned int>(diff.x * 1024);
        y = static_cast<unsigned int>(diff.y * 1024);
        z = static_cast<unsigned int>(diff.z * 1024);
    }

    CU_CALLABLE glm::vec3 from(glm::vec3 const& span, glm::vec3 const& lower) const {
        glm::vec3 pos(x / 1024.f, y / 1024.f, z / 1024.f);
        return (pos * span) + lower;
    }
};

CU_CALLABLE inline float t_compensate(float span) {
    return span / 1024.f / 0.5f;
}

CU_CALLABLE inline box3f leftBounds(box3f const& bounds, float split_pos, float radius, int dim) {
    device::box3f lbounds = bounds;
    lbounds.upper[dim] = split_pos;
    lbounds.upper[dim] += radius + t_compensate(bounds.span()[dim]);
    return lbounds;
}

CU_CALLABLE inline box3f rightBounds(box3f const& bounds, float split_pos, float radius, int dim) {
    device::box3f rbounds = bounds;
    rbounds.lower[dim] = split_pos;
    rbounds.lower[dim] -= radius + t_compensate(bounds.span()[dim]);
    return rbounds;
}
} // namespace device
} // namespace optix_hpg
} // namespace megamol
