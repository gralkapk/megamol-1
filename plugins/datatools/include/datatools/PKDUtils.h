#pragma once

#include <limits>

#include <glm/glm.hpp>

#include "CUDA_ContextCallable.h"

namespace megamol::datatools {
typedef struct box3f {
    box3f()
            : lower{std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                  std::numeric_limits<float>::max()}
            , upper{std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(),
                  std::numeric_limits<float>::lowest()} {}
    ~box3f() = default;
    void extend(glm::vec3 const& pos) {
        lower = glm::vec3(fminf(pos.x, lower.x), fminf(pos.y, lower.y), fminf(pos.z, lower.z));
        upper = glm::vec3(fmaxf(pos.x, upper.x), fmaxf(pos.y, upper.y), fmaxf(pos.z, upper.z));
    }
    void extend(box3f const& box) {
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

struct pkdlet {
    unsigned int begin, end;
    box3f bounds;
};
}
