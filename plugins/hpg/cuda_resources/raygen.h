#pragma once

#include "optix.h"

#include "framestate.h"
#include "perraydata.h"

#include "glm/glm.hpp"

namespace megamol {
namespace hpg {
    namespace optix {
        namespace device {
            struct RayGenData {
                OptixTraversableHandle world;
                int rec_depth;
                glm::uvec2 fbSize;
                glm::vec4* colorBufferPtr;
                FrameState* frameStateBuffer;
            };
        } // namespace device
    }     // namespace optix
} // namespace hpg
} // namespace megamol
