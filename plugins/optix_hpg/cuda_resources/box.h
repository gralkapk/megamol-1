#pragma once

#include <glm/glm.hpp>

#ifndef __CUDACC__
#define CU_CALLABLE
#else
#define CU_CALLABLE __host__ __device__
#endif

namespace megamol {
namespace optix_hpg {
namespace device {

typedef struct box3f {
    CU_CALLABLE box3f()
            : lower{std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                  std::numeric_limits<float>::max()}
            , upper{std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(),
                  std::numeric_limits<float>::lowest()} {}
    CU_CALLABLE ~box3f() = default;
    CU_CALLABLE void extend(glm::vec3 const& pos) {
        lower = glm::vec3(fminf(pos.x, lower.x), fminf(pos.y, lower.y), fminf(pos.z, lower.z));
        upper = glm::vec3(fmaxf(pos.x, upper.x), fmaxf(pos.y, upper.y), fmaxf(pos.z, upper.z));
    }
    CU_CALLABLE void extend(box3f const& box) {
        extend(box.lower);
        extend(box.upper);
    }
    CU_CALLABLE glm::vec3 center() const {
        return (lower + upper) * 0.5f;
    }
    CU_CALLABLE glm::vec3 span() const {
        return upper - lower;
    }
    glm::vec3 lower;
    glm::vec3 upper;
} box3f;

typedef struct box3u32 {
    CU_CALLABLE box3u32()
            : lower{std::numeric_limits<uint32_t>::max(), std::numeric_limits<uint32_t>::max(),
                  std::numeric_limits<uint32_t>::max()}
            , upper{std::numeric_limits<uint32_t>::lowest(), std::numeric_limits<uint32_t>::lowest(),
                  std::numeric_limits<uint32_t>::lowest()} {}
    CU_CALLABLE ~box3u32() = default;
    CU_CALLABLE void extend(glm::uvec3 const& pos) {
        lower = glm::uvec3(fminf(pos.x, lower.x), fminf(pos.y, lower.y), fminf(pos.z, lower.z));
        upper = glm::uvec3(fmaxf(pos.x, upper.x), fmaxf(pos.y, upper.y), fmaxf(pos.z, upper.z));
    }
    CU_CALLABLE void extend(box3u32 const& box) {
        extend(box.lower);
        extend(box.upper);
    }
    CU_CALLABLE glm::uvec3 center() const {
        return glm::uvec3(
            (lower.x >> 2) + (upper.x >> 2), (lower.y >> 2) + (upper.y >> 2), (lower.z >> 2) + (upper.z >> 2));
    }
    CU_CALLABLE glm::uvec3 span() const {
        return upper - lower;
    }
    glm::uvec3 lower;
    glm::uvec3 upper;
} box3u32;
} // namespace device
} // namespace optix_hpg
} // namespace megamol
