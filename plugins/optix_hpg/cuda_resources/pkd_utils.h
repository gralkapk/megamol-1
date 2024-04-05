#pragma once

#include <glm/glm.hpp>

#include "particle.h"

#include "FixedPoint.h"

#ifndef __CUDACC__
#define CU_CALLABLE
#else
#define CU_CALLABLE __host__ __device__
#endif

#define PKD_BOUNDS_CENTER bounds.lower

namespace megamol {
namespace optix_hpg {
inline device::QPKDParticle encode_coord(glm::vec3 const& pos, glm::vec3 const& center, glm::vec3 const& span) {
    #if 0
    constexpr unsigned int digits = 32767u;
    device::QPKDParticle p;
    auto dir = pos - center;
    /*p.sx = dir.x < 0.f;
    p.sy = dir.y < 0.f;
    p.sz = dir.z < 0.f;
    dir = glm::abs(dir);*/
    auto const diff = span / static_cast<float>(digits);
    p.x = static_cast<unsigned int>(dir.x / diff.x);
    p.y = static_cast<unsigned int>(dir.y / diff.y);
    p.z = static_cast<unsigned int>(dir.z / diff.z);

    return p;
    #endif
    decvec3 dec_pos = pos;
    decvec3 dec_center = center;
    auto const diff = dec_pos;
    //-dec_center;
    device::QPKDParticle p;
    p.x = diff.x;
    p.y = diff.y;
    p.z = diff.z;
    return p;
}

inline CU_CALLABLE glm::vec3 decode_coord(
    device::QPKDParticle const& coord /*, glm::vec3 const& center, glm::vec3 const& span*/) {
    #if 0
    constexpr unsigned int digits = 32767u;
    auto diff = span / static_cast<float>(digits);
    /*diff.x = coord.sx ? -diff.x : diff.x;
    diff.y = coord.sy ? -diff.y : diff.y;
    diff.z = coord.sz ? -diff.z : diff.z;*/
#ifndef __CUDACC__
    auto pos = glm::vec3(static_cast<float>(coord.x) * diff.x, static_cast<float>(coord.y) * diff.y,
        static_cast<float>(coord.z) * diff.z);
    pos = pos + center;
#else
    auto const pos = glm::vec3(fmaf(static_cast<float>(coord.x), diff.x, center.x),
        fmaf(static_cast<float>(coord.y), diff.y, center.y), fmaf(static_cast<float>(coord.z), diff.z, center.z));
#endif

    return pos;
    #endif

    decvec3 org = decvec3(coord.x, coord.y, coord.z);
    //decvec3 dec_center = center;
    //auto const pos = org;
    ////+dec_center;
    //glm::vec3 res = pos;
    return org;
}

inline void encode_dim(int dim, device::QPKDParticle& coord) {
    #if 0
    if (dim == 0) {
        coord.dim_x = 1;
    } else {
        coord.dim_x = 0;
    }
    if (dim == 1) {
        coord.dim_y = 1;
    } else {
        coord.dim_y = 0;
    }
    if (dim == 2) {
        coord.dim_z = 1;
    } else {
        coord.dim_z = 0;
    }
    #endif
    coord.dim = dim;
}

inline CU_CALLABLE int decode_dim(device::QPKDParticle const& coord) {
    #if 0
    if (coord.dim_y)
        return 1;
    if (coord.dim_z)
        return 2;
    return 0;
    #endif
    return coord.dim;
}

#ifdef __CUDACC__
namespace device {

inline __device__ PKDParticle const& decode_coord(
    QPKDParticle const& coord /*, glm::vec3 const& center, glm::vec3 const& span*/) {
    //constexpr unsigned int digits = 1023u;
    //auto const diff = span / static_cast<float>(digits);
    ///*auto pos = glm::vec3(static_cast<float>(coord.x) * span.x, static_cast<float>(coord.y) * span.y,
    //    static_cast<float>(coord.z) * span.z);
    //pos = pos + center;*/
    //auto const pos = glm::vec3(fmaf(static_cast<float>(coord.x), diff.x, center.x),
    //    fmaf(static_cast<float>(coord.y), diff.y, center.y), fmaf(static_cast<float>(coord.z), diff.z, center.z));
    PKDParticle p;
    p.dim = megamol::optix_hpg::decode_dim(coord);
    p.pos = megamol::optix_hpg::decode_coord(coord /*, center, span*/);
    return p;
}
} // namespace device
#endif


inline CU_CALLABLE glm::vec3 decode_spart(device::SPKDParticle const& part, device::SPKDlet const& treelet) {
    constexpr const float factor = 1.0f / static_cast<float>(1 << dec_val);
    glm::vec3 pos;
    //device::QPKDParticle qp;
    byte_cast bc;
    bc.ui = 0;
    bc.parts.a = part.x;
    bc.parts.b = treelet.sx[part.sx_idx];
#ifdef __CUDACC__
    pos.x = fmaf(static_cast<float>(bc.ui), factor, treelet.lower.x);
#else
    pos.x = static_cast<float>(bc.ui) * factor;
#endif
    //qp.x = bc.ui;
    bc.parts.a = part.y;
    bc.parts.b = treelet.sy[part.sy_idx];
#ifdef __CUDACC__
    pos.y = fmaf(static_cast<float>(bc.ui), factor, treelet.lower.y);
#else
    pos.y = static_cast<float>(bc.ui) * factor;
#endif
    //qp.y = bc.ui;
    bc.parts.a = part.z;
    bc.parts.b = treelet.sz[part.sz_idx];
#ifdef __CUDACC__
    pos.z = fmaf(static_cast<float>(bc.ui), factor, treelet.lower.z);
#else
    pos.z = static_cast<float>(bc.ui) * factor;
#endif
    //qp.z = bc.ui;
    //return megamol::optix_hpg::decode_coord(qp /*, glm::vec3(), glm::vec3()*/);
    return pos;
}

} // namespace optix_hpg
} // namespace megamol
