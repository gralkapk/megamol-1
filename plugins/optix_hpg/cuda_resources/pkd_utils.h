#pragma once

#include <glm/glm.hpp>

#include "particle.h"

#ifndef __CUDACC__
#define CU_CALLABLE
#else
#define CU_CALLABLE __host__ __device__
#endif

namespace megamol {
namespace optix_hpg {
inline device::QPKDParticle encode_coord(glm::vec3 const& pos, glm::vec3 const& center, glm::vec3 const& span) {
    constexpr unsigned int digits = 511u;
    device::QPKDParticle p;
    auto dir = pos - center;
    p.sx = dir.x < 0.f;
    p.sy = dir.y < 0.f;
    p.sz = dir.z < 0.f;
    dir = glm::abs(dir);
    auto const diff = span / static_cast<float>(digits);
    p.x = static_cast<unsigned int>(dir.x / diff.x);
    p.y = static_cast<unsigned int>(dir.y / diff.y);
    p.z = static_cast<unsigned int>(dir.z / diff.z);

    return p;
}

inline CU_CALLABLE glm::vec3 decode_coord(device::QPKDParticle const& coord, glm::vec3 const& center, glm::vec3 const& span) {
    constexpr unsigned int digits = 511u;
    auto diff = span / static_cast<float>(digits);
    diff.x = coord.sx ? -diff.x : diff.x;
    diff.y = coord.sy ? -diff.y : diff.y;
    diff.z = coord.sz ? -diff.z : diff.z;
#ifndef __CUDACC__
    auto pos = glm::vec3(static_cast<float>(coord.x) * diff.x, static_cast<float>(coord.y) * diff.y,
        static_cast<float>(coord.z) * diff.z);
    pos = pos + center;
#else
    auto const pos = glm::vec3(fmaf(static_cast<float>(coord.x), diff.x, center.x),
        fmaf(static_cast<float>(coord.y), diff.y, center.y), fmaf(static_cast<float>(coord.z), diff.z, center.z));
#endif

    return pos;
}

#ifdef __CUDACC__
namespace device {

inline __device__ PKDParticle const& decode_coord(
    QPKDParticle const& coord, glm::vec3 const& center, glm::vec3 const& span) {
    //constexpr unsigned int digits = 1023u;
    //auto const diff = span / static_cast<float>(digits);
    ///*auto pos = glm::vec3(static_cast<float>(coord.x) * span.x, static_cast<float>(coord.y) * span.y,
    //    static_cast<float>(coord.z) * span.z);
    //pos = pos + center;*/
    //auto const pos = glm::vec3(fmaf(static_cast<float>(coord.x), diff.x, center.x),
    //    fmaf(static_cast<float>(coord.y), diff.y, center.y), fmaf(static_cast<float>(coord.z), diff.z, center.z));
    PKDParticle p;
    p.dim = coord.dim;
    p.pos = megamol::optix_hpg::decode_coord(coord, center, span);
    return p;
}
} // namespace device
#endif
} // namespace optix_hpg
} // namespace megamol
