#pragma once

#include <glm/glm.hpp>

namespace megamol {
namespace optix_hpg {
namespace device {
struct Particle {
    //float x, y, z;
    glm::vec3 pos;
};

struct PKDParticle {
    glm::vec3 pos;
    float dim;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
