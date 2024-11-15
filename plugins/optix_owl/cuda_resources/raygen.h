#pragma once

#include <owl/common/math/vec.h>

#include "framestate.h"
#include "particle.h"

namespace megamol {
namespace optix_owl {
namespace device {
using namespace owl::common;
struct RayGenData {
    OptixTraversableHandle world;
    int rec_depth;
    vec2ui fbSize;
    uint32_t* colorBufferPtr;
    float4* accumBufferPtr;
    FrameState* frameStateBuffer;
    float intensity;
    vec3f background;
};
} // namespace device
} // namespace optix_owl
} // namespace megamol
