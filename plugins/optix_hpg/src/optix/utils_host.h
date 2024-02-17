#pragma once

#include <glm/glm.hpp>

#define MM_OPTIX_RAYGEN_ANNOTATION __raygen__

#define MM_OPTIX_INTERSECTION_ANNOTATION __intersection__

#define MM_OPTIX_ANYHIT_ANNOTATION __anyhit__

#define MM_OPTIX_CLOSESTHIT_ANNOTATION __closesthit__

#define MM_OPTIX_MISS_ANNOTATION __miss__

#define MM_OPTIX_DIRECT_CALLABLE_ANNOTATION __direct_callable__

#define MM_OPTIX_CONTINUATION_CALLABLE_ANNOTATION __continuation_callable__

#define MM_OPTIX_EXCEPTION_ANNOTATION __exception__

#define MM_OPTIX_BOUNDS_ANNOTATION __boundsKernel__


#define MM_OPTIX_RAYGEN_ANNOTATION_STRING "__raygen__"

#define MM_OPTIX_INTERSECTION_ANNOTATION_STRING "__intersection__"

#define MM_OPTIX_ANYHIT_ANNOTATION_STRING "__anyhit__"

#define MM_OPTIX_CLOSESTHIT_ANNOTATION_STRING "__closesthit__"

#define MM_OPTIX_MISS_ANNOTATION_STRING "__miss__"

#define MM_OPTIX_DIRECT_CALLABLE_ANNOTATION_STRING "__direct_callable__"

#define MM_OPTIX_CONTINUATION_CALLABLE_ANNOTATION_STRING "__continuation_callable__"

#define MM_OPTIX_EXCEPTION_ANNOTATION_STRING "__exception__"

#define MM_OPTIX_BOUNDS_ANNOTATION_STRING "__boundsKernel__"

namespace megamol {
namespace optix_hpg {

typedef struct box3f {
    box3f()
            : lower{std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                  std::numeric_limits<float>::max()}
            , upper{std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(),
                  std::numeric_limits<float>::lowest()} {}
    void extend(glm::vec3 const& pos) {
        lower = glm::vec3(fminf(pos.x, lower.x), fminf(pos.y, lower.y), fminf(pos.z, lower.z));
        upper = glm::vec3(fmaxf(pos.x, upper.x), fmaxf(pos.y, upper.y), fmaxf(pos.z, upper.z));
    }
    void extend(box3f const& box) {
        extend(box.lower);
        extend(box.upper);
    }
    glm::vec3 center() const {
        return (lower + upper) * 0.5f;
    }
    glm::vec3 span() const {
        return upper - lower;
    }
    glm::vec3 lower;
    glm::vec3 upper;
} box3f;

typedef struct RayH {
    RayH(glm::vec3 const& org, glm::vec3 const& dir, float tmin, float tmax)
            : origin(org)
            , direction(dir)
            , tMin(tmin)
            , tMax(tmax) {}
    glm::vec3 origin;
    glm::vec3 direction;
    float tMin;
    float tMax;
};
} // namespace optix_hpg
} // namespace megamol
