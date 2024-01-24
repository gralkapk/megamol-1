#pragma once

#include "particle.h"
#include "perraydata.h"

namespace megamol {
namespace optix_hpg {
namespace device {
struct PKDGeoData {
    PKDParticle* particleBufferPtr;
    float* radiusBufferPtr;
    glm::vec4* colorBufferPtr;
    bool hasGlobalRadius;
    float radius;
    bool hasColorData;
    glm::vec4 globalColor;
    unsigned int particleCount;
    glm::vec3 worldBounds_lower;
    glm::vec3 worldBounds_upper;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
