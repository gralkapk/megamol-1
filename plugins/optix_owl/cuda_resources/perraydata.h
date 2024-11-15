#pragma once

#include <owl/common/math/vec.h>

namespace megamol {
namespace optix_owl {
namespace device {
struct PerRayData {
    int particleID;
    float t;
    owl::vec3f pos;
    owl::vec3f color;
};
} // namespace device
} // namespace optix_owl
} // namespace megamol
