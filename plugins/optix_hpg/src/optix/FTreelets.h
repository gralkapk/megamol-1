#pragma once

#include <vector>

#include "particle.h"

namespace megamol::optix_hpg {
inline device::FPKDParticle convert_to_fparticle(
    device::PKDParticle const& particle, glm::dvec3 const& span, glm::dvec3 const& lower) {
    constexpr uint64_t const factor = 1 << 23;
    auto const org_pos = glm::dvec3(particle.pos);
    auto new_pos = (org_pos - lower) / span;
    new_pos *= factor;
    device::FPKDParticle ret;
    ret.pos = new_pos;
    return ret;
}

inline std::vector<device::FPKDParticle> convert_to_fparticles(
    std::vector<device::PKDParticle> const& particles, device::box3f const& bounds) {
    std::vector<device::FPKDParticle> ret(particles.size());

    auto const span = glm::dvec3(bounds.span());
    auto const lower = glm::dvec3(bounds.lower);

    for (size_t i = 0; i < particles.size(); ++i) {
        ret[i] = convert_to_fparticle(particles[i], span, lower);
    }

    return ret;
}
} // namespace megamol::optix_hpg
