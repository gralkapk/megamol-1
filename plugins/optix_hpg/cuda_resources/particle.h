#pragma once

#include <glm/glm.hpp>

#include "optix/utils_host.h"

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

struct PKDlet {
    //! bounding box of all particles (including the radius)
    box3f bounds;
    //! begin/end range in the common particles array
    size_t begin, end;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
