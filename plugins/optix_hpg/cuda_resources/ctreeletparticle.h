#pragma once

#include "box.h"

namespace megamol {
namespace optix_hpg {
namespace device {
struct CPKDlet {
    box3f bounds;
    box3u32 fbounds;
    unsigned int begin, end;
    char offset[3];
    unsigned short prefix[3];
};

struct CPKDParticle {
    unsigned int dim : 2;
    unsigned int x : 10;
    unsigned int y : 10;
    unsigned int z : 10;

    CU_CALLABLE glm::vec3 from(CPKDlet const& tl, glm::vec3 const& span, glm::vec3 const& lower) const {
        glm::vec3 pos(x + tl.fbounds.lower.x, y + tl.fbounds.lower.y, z + tl.fbounds.lower.z);

        pos.x += static_cast<uint32_t>(tl.prefix[0]) << tl.offset[0];
        pos.y += static_cast<uint32_t>(tl.prefix[1]) << tl.offset[1];
        pos.z += static_cast<uint32_t>(tl.prefix[2]) << tl.offset[2];

        constexpr uint64_t const factor = 1 << 23;
        pos.x /= static_cast<float>(factor);
        pos.y /= static_cast<float>(factor);
        pos.z /= static_cast<float>(factor);

        pos = pos * span + lower;
        //pos = pos * span + lower;

        return pos;
    }
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
