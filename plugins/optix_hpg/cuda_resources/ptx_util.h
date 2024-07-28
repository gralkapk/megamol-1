#pragma once

#include <glm/glm.hpp>

namespace megamol {
namespace optix_hpg {
namespace device {

inline __device__ bool clipToBounds(const Ray& ray, const box3f& bounds, float& t0, float& t1) {
    glm::vec3 t_lower = (bounds.lower - ray.origin) / ray.direction;
    glm::vec3 t_upper = (bounds.upper - ray.origin) / ray.direction;

    glm::vec3 t_min3 = min(t_lower, t_upper);
    glm::vec3 t_max3 = max(t_lower, t_upper);

    t0 = fmaxf(ray.tmin, glm::max(t_min3.x, glm::max(t_min3.y, t_min3.z)));
    t1 = fminf(ray.tmax, glm::min(t_max3.x, glm::min(t_max3.y, t_max3.z)));
    return t0 < t1;
}

inline __device__ bool clipToBounds(const Ray& ray, const datatools::box3f& bounds, float& t0, float& t1) {
    glm::vec3 t_lower = (bounds.lower - ray.origin) / ray.direction;
    glm::vec3 t_upper = (bounds.upper - ray.origin) / ray.direction;

    glm::vec3 t_min3 = min(t_lower, t_upper);
    glm::vec3 t_max3 = max(t_lower, t_upper);

    t0 = fmaxf(ray.tmin, glm::max(t_min3.x, glm::max(t_min3.y, t_min3.z)));
    t1 = fminf(ray.tmax, glm::min(t_max3.x, glm::min(t_max3.y, t_max3.z)));
    return t0 < t1;
}

inline __device__ int getDim(glm::vec3 const& pos) {
    return *((unsigned int*) &pos.x) & 3;
}

inline __device__ bool intersectSphere(
    const PKDParticle& particle, const float particleRadius, const Ray& ray, float& hit_t) {
    // Raytracing Gems Intersection Code (Chapter 7)
    //const glm::vec3 pos = glm::vec3(particle.x, particle.y, particle.z);
    const glm::vec3 oc = ray.origin - particle.pos;
    const float sqrRad = particleRadius * particleRadius;

    // const float  a = dot(ray.direction, ray.direction);
    const float b = glm::dot(-oc, ray.direction);
    const glm::vec3 temp = oc + b * ray.direction;
    const float delta = sqrRad - glm::dot(temp, temp);

    if (delta < 0.0f)
        return false;

    const float c = glm::dot(oc, oc) - sqrRad;
    const float q = b + copysignf(sqrtf(delta), b);

    {
        float temp = fminf(c / q, q);
        if (temp < hit_t && temp > ray.tmin) {
            hit_t = temp;
            return true;
        }
    }

    return false;
}

inline __device__ bool intersectSphere(const glm::vec3& pos, const float particleRadius, const Ray& ray, float& hit_t) {
    // Raytracing Gems Intersection Code (Chapter 7)
    //const glm::vec3 pos = glm::vec3(particle.x, particle.y, particle.z);
    const glm::vec3 oc = ray.origin - pos;
    const float sqrRad = particleRadius * particleRadius;

    // const float  a = dot(ray.direction, ray.direction);
    const float b = glm::dot(-oc, ray.direction);
    const glm::vec3 temp = oc + b * ray.direction;
    const float delta = sqrRad - glm::dot(temp, temp);

    if (delta < 0.0f)
        return false;

    const float c = glm::dot(oc, oc) - sqrRad;
    const float q = b + copysignf(sqrtf(delta), b);

    {
        float temp = fminf(c / q, q);
        if (temp < hit_t && temp > ray.tmin) {
            hit_t = temp;
            return true;
        }
    }

    return false;
}

struct StackEntry {
    float t0, t1;
    unsigned int nodeID;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
