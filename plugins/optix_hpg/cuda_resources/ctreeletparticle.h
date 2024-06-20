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
    morton_prefix_t prefix;
};

struct C2PKDParticle {
    unsigned int dim : 2;
    unsigned int code : 30;

    CU_CALLABLE glm::vec3 from(morton_prefix_t const prefix, glm::vec3 const& span, glm::vec3 const& lower) const {
        auto const combined_code = (static_cast<uint64_t>(code) << MortonConfig::code_offset) +
                                   (static_cast<uint64_t>(prefix) << MortonConfig::prefix_offset);

        static float const ffactor = MortonConfig::factor;

        uint32_t x, y, z;
        morton_decode(combined_code, x, y, z);
#ifdef __CUDACC__
        glm::vec3 basePos(
            fmaf(x / ffactor, span.x, lower.x), fmaf(y / ffactor, span.y, lower.y), fmaf(z / ffactor, span.z, lower.z));
#else
        glm::vec3 basePos(x / ffactor, y / ffactor, z / ffactor);
        basePos *= span;
        basePos += lower;
#endif

        return basePos;
    }
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
