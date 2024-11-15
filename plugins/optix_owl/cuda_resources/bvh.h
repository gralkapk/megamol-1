#pragma once

#include <owl/common/math/box.h>
#include <owl/common/math/vec.h>

#include "particle.h"

namespace megamol {
namespace optix_owl {
namespace device {
using namespace owl::common;
struct BVHGeomData {
    Particle* particleBuffer;
    vec3f* colorBuffer;
    float particleRadius;
};
} // namespace device
} // namespace optix_owl
} // namespace megamol
