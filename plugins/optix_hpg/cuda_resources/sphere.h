#pragma once

#include "particle.h"
#include "perraydata.h"

namespace megamol {
namespace optix_hpg {
namespace device {
struct SphereGeoData {
    Particle* particleBufferPtr;
    float* radiusBufferPtr;
    color_t* colorBufferPtr;
    bool hasGlobalRadius;
    float radius;
    bool hasColorData;
    color_t globalColor;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
