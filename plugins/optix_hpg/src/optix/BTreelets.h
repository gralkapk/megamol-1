#pragma once

#include <cstdint>

#include "PKDUtils.h"
#include "box.h"
#include "particle.h"

namespace megamol::optix_hpg {
inline void convert_blets(uint64_t P, uint64_t N, device::PKDParticle const* particles,
    device::BTParticle* out_particles, float radius, device::box3f bounds) {
    if (P >= N)
        return;

    auto& out_par = out_particles[P];
    out_par.dim = particles[P].dim;

    auto const span = bounds.span();

    out_par.to(particles[P].pos, span, bounds.lower);
    /*auto const diff = (particles[P].pos - bounds.lower) / span;
    out_par.x = static_cast<unsigned int>(diff.x * 1024);
    out_par.y = static_cast<unsigned int>(diff.y * 1024);
    out_par.z = static_cast<unsigned int>(diff.z * 1024);*/

    auto const split_pos = (out_par.get(out_par.dim) / 1023.f) * span[out_par.dim] + bounds.lower[out_par.dim];
    //auto const split_pos = particles[P].pos[out_par.dim];

    auto const lbounds = device::leftBounds(bounds, split_pos, radius, out_par.dim);
    auto const rbounds = device::rightBounds(bounds, split_pos, radius, out_par.dim);

    /*device::box3f lbounds = bounds;
    lbounds.upper[out_par.dim] = split_pos;
    lbounds.upper[out_par.dim] += radius + device::t_compensate(span[out_par.dim]);

    device::box3f rbounds = bounds;
    rbounds.lower[out_par.dim] = split_pos;
    rbounds.lower[out_par.dim] -= radius + device::t_compensate(span[out_par.dim]);*/

    convert_blets(lChild(P), N, particles, out_particles, radius, lbounds);
    convert_blets(rChild(P), N, particles, out_particles, radius, rbounds);
}

inline void reconstruct_blets(uint64_t P, uint64_t N, device::PKDParticle const* particles,
    device::BTParticle* out_particles, float radius, device::box3f bounds, glm::vec3* original, glm::vec3* rec,
    glm::vec3* diffs) {
    if (P >= N)
        return;

    auto& out_par = out_particles[P];

    auto const span = bounds.span();

    auto const pos = out_par.from(span, bounds.lower);

    auto const split_pos = pos[out_par.dim];

    auto const lbounds = device::leftBounds(bounds, split_pos, radius, out_par.dim);
    auto const rbounds = device::rightBounds(bounds, split_pos, radius, out_par.dim);

    glm::dvec3 org_pos = particles[P].pos;
    glm::dvec3 new_pos = pos;
    auto const diff = new_pos - org_pos;

    original[P] = org_pos;
    rec[P] = new_pos;
    diffs[P] = diff;

    reconstruct_blets(lChild(P), N, particles, out_particles, radius, lbounds, original, rec, diffs);
    reconstruct_blets(rChild(P), N, particles, out_particles, radius, rbounds, original, rec, diffs);
}
} // namespace megamol::optix_hpg