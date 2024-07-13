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
    color_t* colorBufferPtr;
    bool hasGlobalRadius;
    float radius;
    bool hasColorData;
    color_t globalColor;
    unsigned int particleCount;
    box3f worldBounds;
};
struct TreeletsGeoData {
    PKDParticle* particleBufferPtr;
    float* radiusBufferPtr;
    color_t* colorBufferPtr;
    PKDlet* treeletBufferPtr;
    bool hasGlobalRadius;
    float radius;
    bool hasColorData;
    color_t globalColor;
    unsigned int particleCount;
    box3f worldBounds;
};
struct QTreeletsGeoData {
    QPKDParticle* particleBufferPtr;
    color_t* colorBufferPtr;
    PKDlet* treeletBufferPtr;
    float radius;
    bool hasColorData;
    color_t globalColor;
    unsigned int particleCount;
};
struct STreeletsGeoData {
    SPKDParticle* particleBufferPtr;
    color_t* colorBufferPtr;
    SPKDlet* treeletBufferPtr;
    float radius;
    bool hasColorData;
    color_t globalColor;
    unsigned int particleCount;
};
struct QPKDTreeletsGeoData {
    void* particleBufferPtr;
    color_t* colorBufferPtr;
    QPKDlet* treeletBufferPtr;
    char* expXBuffer;
    char* expYBuffer;
    char* expZBuffer;
    float radius;
    bool hasColorData;
    color_t globalColor;
    unsigned int particleCount;
    int selectedType;
    char use_localtables;
};
struct BTreeletsGeoData {
    BTParticle* particleBufferPtr;
    color_t* colorBufferPtr;
    PKDlet* treeletBufferPtr;
    float radius;
    bool hasColorData;
    color_t globalColor;
    unsigned int particleCount;
};
struct CTreeletsGeoData {
    box3f bounds;
    C2PKDParticle* particleBufferPtr;
    color_t* colorBufferPtr;
    C2PKDlet* treeletBufferPtr;
    float radius;
    bool hasColorData;
    color_t globalColor;
    unsigned int particleCount;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
