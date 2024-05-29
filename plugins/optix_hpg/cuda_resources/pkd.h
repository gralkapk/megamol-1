#pragma once

#include "box.h"
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
struct QTreeletsGeoData {
    QPKDParticle* particleBufferPtr;
    glm::vec4* colorBufferPtr;
    PKDlet* treeletBufferPtr;
    float radius;
    bool hasColorData;
    glm::vec4 globalColor;
    unsigned int particleCount;
};
struct STreeletsGeoData {
    SPKDParticle* particleBufferPtr;
    glm::vec4* colorBufferPtr;
    SPKDlet* treeletBufferPtr;
    float radius;
    bool hasColorData;
    glm::vec4 globalColor;
    unsigned int particleCount;
};
struct QPKDTreeletsGeoData {
    void* particleBufferPtr;
    glm::vec4* colorBufferPtr;
    QPKDlet* treeletBufferPtr;
    char* expXBuffer;
    char* expYBuffer;
    char* expZBuffer;
    float radius;
    bool hasColorData;
    glm::vec4 globalColor;
    unsigned int particleCount;
    int selectedType;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
