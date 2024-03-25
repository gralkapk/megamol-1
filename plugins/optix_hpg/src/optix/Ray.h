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
typedef struct Ray {
    CU_CALLABLE Ray(float3 const& org, float3 const& dir, float tmin, float tmax)
            : origin(org.x, org.y, org.z)
            , direction(dir.x, dir.y, dir.z)
            , tmin(tmin)
            , tmax(tmax) {}
    CU_CALLABLE Ray(glm::vec3 const& org, glm::vec3 const& dir, float tmin, float tmax)
            : origin(org)
            , direction(dir)
            , tmin(tmin)
            , tmax(tmax) {}

    glm::vec3 origin;
    glm::vec3 direction;
    float tmin;
    float tmax;
} Ray;
} // namespace device
} // namespace optix_hpg
} // namespace megamol
