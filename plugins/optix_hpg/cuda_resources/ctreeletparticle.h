#pragma once

#include "box.h"
#include "morton_util.h"

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

struct C2PKDlet {
    box3f bounds;
    unsigned int begin, end;
    unsigned short prefix;
};

struct C2PKDParticle {
    unsigned int dim : 2;
    unsigned int code : 30;

    CU_CALLABLE glm::vec3 from(unsigned short const prefix, glm::vec3 const& span, glm::vec3 const& lower,
        int const code_offset, int const prefix_offset, float const factor) const {
        auto const combined_code =
            (static_cast<uint64_t>(code) << code_offset) + (static_cast<uint64_t>(prefix) << prefix_offset);

        uint32_t x, y, z;
        morton_decode(combined_code, x, y, z);
        glm::vec3 basePos(x / factor, y / factor, z / factor);
        basePos *= span;
        basePos += lower;

        return basePos;
    }
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
