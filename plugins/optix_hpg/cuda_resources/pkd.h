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
    box3f worldBounds;
};
struct TreeletsGeoData {
    PKDParticle* particleBufferPtr;
    float* radiusBufferPtr;
    glm::vec4* colorBufferPtr;
    PKDlet* treeletBufferPtr;
    bool hasGlobalRadius;
    float radius;
    bool hasColorData;
    glm::vec4 globalColor;
    unsigned int particleCount;
    box3f worldBounds;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
