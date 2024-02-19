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
} // namespace device
} // namespace optix_hpg
} // namespace megamol
