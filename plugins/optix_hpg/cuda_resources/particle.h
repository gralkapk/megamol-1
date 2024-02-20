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

struct QPKDParticle {
    #if 0
    unsigned int dim : 2;
    unsigned int sx : 1;
    unsigned int x : 9;
    unsigned int sy : 1;
    unsigned int y : 9;
    unsigned int sz : 1;
    unsigned int z : 9;
    #endif
    unsigned short dim_x : 1;
    unsigned short x : 15;
    unsigned short dim_y : 1;
    unsigned short y : 15;
    unsigned short dim_z : 1;
    unsigned short z : 15;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
