// PKD implementation:
// ======================================================================== //
// Copyright 2018-2019 Ingo Wald                                            //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

// ======================================================================== //
// Modified 2019-2020 VISUS - University of Stuttgart                       //
// ======================================================================== //

#include "pkd.h"

#include "optix/utils_device.h"
#include "pkd_utils.h"

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

struct StackEntry {
    float t0, t1;
    unsigned int nodeID;
};

MM_OPTIX_INTERSECTION_KERNEL(pkd_intersect)() {
    const auto& self = getProgramData<PKDGeoData>();

    float t0, t1;
    {
        auto const ray =
            Ray(optixGetWorldRayOrigin(), optixGetWorldRayDirection(), optixGetRayTmin(), optixGetRayTmax());
        if (!clipToBounds(ray, self.worldBounds, t0, t1))
            return;
    }
    auto const ray = Ray(optixGetWorldRayOrigin(), optixGetWorldRayDirection(), fmaxf(optixGetRayTmin(), t0),
        fminf(optixGetRayTmax(), t1));

    int nodeID = 0;
    float tmp_hit_t = t1;
    int tmp_hit_primID = -1;

    enum { STACK_DEPTH = 32 };
    StackEntry stackBase[STACK_DEPTH];
    StackEntry* stackPtr = stackBase;

    const int dir_sign[3] = {ray.direction.x < 0.f, ray.direction.y < 0.f, ray.direction.z < 0.f};
    const float org[3] = {ray.origin.x, ray.origin.y, ray.origin.z};
    const float rdir[3] = {
        (fabsf(ray.direction.x) <= 1e-8f) ? 1e8f : 1.f / ray.direction.x,
        (fabsf(ray.direction.y) <= 1e-8f) ? 1e8f : 1.f / ray.direction.y,
        (fabsf(ray.direction.z) <= 1e-8f) ? 1e8f : 1.f / ray.direction.z,
    };
    unsigned int const numParticles = self.particleCount;
    float const particleRadius = self.radius;

    while (1) {
        // while we have anything to traverse ...

        while (1) {
            // while we can go down
            const PKDParticle& particle = self.particleBufferPtr[nodeID];
            int const dim = particle.dim;

            const float t_slab_lo = (particle.pos[dim] - particleRadius - org[dim]) * rdir[dim];
            const float t_slab_hi = (particle.pos[dim] + particleRadius - org[dim]) * rdir[dim];

            const float t_slab_nr = fminf(t_slab_lo, t_slab_hi);
            const float t_slab_fr = fmaxf(t_slab_lo, t_slab_hi);

            // -------------------------------------------------------
            // compute potential sphere interval, and intersect if necessary
            // -------------------------------------------------------
            const float sphere_t0 = fmaxf(t0, t_slab_nr);
            const float sphere_t1 = fminf(fminf(t_slab_fr, t1), tmp_hit_t);

            if (sphere_t0 < sphere_t1) {
                if (intersectSphere(particle, particleRadius, ray, tmp_hit_t))
                    tmp_hit_primID = nodeID;
            }

            // -------------------------------------------------------
            // compute near and far side intervals
            // -------------------------------------------------------
            const float nearSide_t0 = t0;
            const float nearSide_t1 = fminf(fminf(t_slab_fr, t1), tmp_hit_t);

            const float farSide_t0 = fmaxf(t0, t_slab_nr);
            const float farSide_t1 = fminf(t1, tmp_hit_t);


            // -------------------------------------------------------
            // logic
            // -------------------------------------------------------
            const int nearSide_nodeID = 2 * nodeID + 1 + dir_sign[dim];
            const int farSide_nodeID = 2 * nodeID + 2 - dir_sign[dim];

            const bool nearSide_valid = nearSide_nodeID < numParticles;
            const bool farSide_valid = farSide_nodeID < numParticles;

            const bool need_nearSide = nearSide_valid && nearSide_t0 < nearSide_t1;
            const bool need_farSide = farSide_valid && farSide_t0 < farSide_t1;

            if (!(need_nearSide || need_farSide))
                break; // pop ...

            if (need_nearSide && need_farSide) {
                stackPtr->t0 = farSide_t0;
                stackPtr->t1 = farSide_t1;
                stackPtr->nodeID = farSide_nodeID;
                ++stackPtr;

                nodeID = nearSide_nodeID;
                t0 = nearSide_t0;
                t1 = nearSide_t1;
                continue;
            }

            nodeID = need_nearSide ? nearSide_nodeID : farSide_nodeID;
            t0 = need_nearSide ? nearSide_t0 : farSide_t0;
            t1 = need_nearSide ? nearSide_t1 : farSide_t1;
        }
        // -------------------------------------------------------
        // pop
        // -------------------------------------------------------
        while (1) {
            if (stackPtr == stackBase) {
                // can't pop any more - done.
                if (tmp_hit_primID >= 0 && tmp_hit_t < ray.tmax) {
                    optixReportIntersection(tmp_hit_t, 0, tmp_hit_primID);
                }
                return;
            }
            --stackPtr;
            t0 = stackPtr->t0;
            t1 = fminf(stackPtr->t1, tmp_hit_t);
            nodeID = stackPtr->nodeID;
            if (t1 <= t0)
                continue;
            break;
        }
    }
}


MM_OPTIX_CLOSESTHIT_KERNEL(pkd_closesthit)() {
    //const int primID = optixGetPrimitiveIndex();
    const unsigned int primID = optixGetAttribute_0();
    PerRayData& prd = getPerRayData<PerRayData>();

    const auto& self = getProgramData<PKDGeoData>();

    prd.particleID = primID;
    const PKDParticle& particle = self.particleBufferPtr[primID];
    prd.pos = particle.pos;
    glm::vec3 geo_col = glm::vec3(self.globalColor);
    if (self.hasColorData) {
        geo_col = glm::vec3(self.colorBufferPtr[primID]);
    }
    prd.albedo = geo_col;
    prd.t = optixGetRayTmax();
    set_depth(prd, optixGetRayTmax());
}


MM_OPTIX_CLOSESTHIT_KERNEL(pkd_closesthit_occlusion)() {
    optixSetPayload_0(1);
}


MM_OPTIX_BOUNDS_KERNEL(pkd_bounds)
(const void* geomData, const float* radData, float radius, box3f& primBounds, const unsigned int primID) {}


MM_OPTIX_INTERSECTION_KERNEL(treelets_intersect)
() {
    const int treeletID = optixGetPrimitiveIndex();
    const auto& self = getProgramData<TreeletsGeoData>();
    const auto treelet = self.treeletBufferPtr[treeletID];

    const int begin = treelet.begin;
    const int size = treelet.end - begin;
    {
        float t0, t1;
        {
            auto const ray =
                Ray(optixGetWorldRayOrigin(), optixGetWorldRayDirection(), optixGetRayTmin(), optixGetRayTmax());
            if (!clipToBounds(ray, treelet.bounds, t0, t1))
                return;
        }
        auto const ray = Ray(optixGetWorldRayOrigin(), optixGetWorldRayDirection(), fmaxf(optixGetRayTmin(), t0),
            fminf(optixGetRayTmax(), t1));

        int nodeID = 0;
        float tmp_hit_t = ray.tmax;
        int tmp_hit_primID = -1;

        enum { STACK_DEPTH = 12 };
        StackEntry stackBase[STACK_DEPTH];
        StackEntry* stackPtr = stackBase;

        const int dir_sign[3] = {ray.direction.x < 0.f, ray.direction.y < 0.f, ray.direction.z < 0.f};
        const float org[3] = {ray.origin.x, ray.origin.y, ray.origin.z};
        const float rdir[3] = {
            (fabsf(ray.direction.x) <= 1e-8f) ? 1e8f : 1.f / ray.direction.x,
            (fabsf(ray.direction.y) <= 1e-8f) ? 1e8f : 1.f / ray.direction.y,
            (fabsf(ray.direction.z) <= 1e-8f) ? 1e8f : 1.f / ray.direction.z,
        };

        while (1) {
            // while we have anything to traverse ...

            while (1) {
                // while we can go down

                const int particleID = nodeID + begin;
                const PKDParticle& particle = self.particleBufferPtr[particleID];
                int const dim = particle.dim;

                const float t_slab_lo = (particle.pos[dim] - self.radius - org[dim]) * rdir[dim];
                const float t_slab_hi = (particle.pos[dim] + self.radius - org[dim]) * rdir[dim];

                const float t_slab_nr = fminf(t_slab_lo, t_slab_hi);
                const float t_slab_fr = fmaxf(t_slab_lo, t_slab_hi);

                // -------------------------------------------------------
                // compute potential sphere interval, and intersect if necessary
                // -------------------------------------------------------
                const float sphere_t0 = fmaxf(t0, t_slab_nr);
                const float sphere_t1 = fminf(fminf(t_slab_fr, t1), tmp_hit_t);

                if (sphere_t0 < sphere_t1) {
                    if (intersectSphere(particle, self.radius, ray, tmp_hit_t)) {
                        tmp_hit_primID = particleID;
                    }
                }

                // -------------------------------------------------------
                // compute near and far side intervals
                // -------------------------------------------------------
                const float nearSide_t0 = t0;
                const float nearSide_t1 = fminf(fminf(t_slab_fr, t1), tmp_hit_t);

                const float farSide_t0 = fmaxf(t0, t_slab_nr);
                const float farSide_t1 = fminf(t1, tmp_hit_t);

                // -------------------------------------------------------
                // logic
                // -------------------------------------------------------
                const int nearSide_nodeID = 2 * nodeID + 1 + dir_sign[dim];
                const int farSide_nodeID = 2 * nodeID + 2 - dir_sign[dim];

                const bool nearSide_valid = nearSide_nodeID < size;
                const bool farSide_valid = farSide_nodeID < size;

                const bool need_nearSide = nearSide_valid && nearSide_t0 < nearSide_t1;
                const bool need_farSide = farSide_valid && farSide_t0 < farSide_t1;

                if (!(need_nearSide || need_farSide))
                    break; // pop ...

                if (need_nearSide && need_farSide) {
                    stackPtr->nodeID = farSide_nodeID;
                    stackPtr->t0 = farSide_t0;
                    stackPtr->t1 = farSide_t1;
                    ++stackPtr;

                    nodeID = nearSide_nodeID;
                    t0 = nearSide_t0;
                    t1 = nearSide_t1;
                    continue;
                }

                nodeID = need_nearSide ? nearSide_nodeID : farSide_nodeID;
                t0 = need_nearSide ? nearSide_t0 : farSide_t0;
                t1 = need_nearSide ? nearSide_t1 : farSide_t1;
            }
            // -------------------------------------------------------
            // pop
            // -------------------------------------------------------
            while (1) {
                if (stackPtr == stackBase) {
                    // can't pop any more - done.
                    if (tmp_hit_primID >= 0 && tmp_hit_t < ray.tmax) {
                        optixReportIntersection(tmp_hit_t, 0, tmp_hit_primID);
                    }
                    return;
                }
                --stackPtr;
                t0 = stackPtr->t0;
                t1 = stackPtr->t1;
                nodeID = stackPtr->nodeID;
                t1 = fminf(t1, tmp_hit_t);
                if (t1 <= t0)
                    continue;
                break;
            }
        }
    }
}


MM_OPTIX_CLOSESTHIT_KERNEL(treelets_closesthit)
() {
    const unsigned int primID = optixGetAttribute_0();
    PerRayData& prd = getPerRayData<PerRayData>();

    const auto& self = getProgramData<TreeletsGeoData>();

    prd.particleID = primID;
    const PKDParticle& particle = self.particleBufferPtr[primID];
    prd.pos = particle.pos;
    glm::vec3 geo_col = glm::vec3(self.globalColor);
    if (self.hasColorData) {
        geo_col = glm::vec3(self.colorBufferPtr[primID]);
    }
    prd.albedo = geo_col;
    prd.t = optixGetRayTmax();
    set_depth(prd, optixGetRayTmax());
}


MM_OPTIX_CLOSESTHIT_KERNEL(treelets_closesthit_occlusion)
() {
    optixSetPayload_0(1);
}


MM_OPTIX_BOUNDS_KERNEL(treelets_bounds)
(const void* geomData, const float* radData, float radius, box3f& primBounds, const unsigned int primID) {}


//inline __device__ PKDParticle const& decode_coord(QPKDParticle const& coord, glm::vec3 const& center, glm::vec3 const& span) {
//    constexpr unsigned int digits = 1023u;
//    auto const diff = span / static_cast<float>(digits);
//    /*auto pos = glm::vec3(static_cast<float>(coord.x) * span.x, static_cast<float>(coord.y) * span.y,
//        static_cast<float>(coord.z) * span.z);
//    pos = pos + center;*/
//    auto const pos = glm::vec3(fmaf(static_cast<float>(coord.x), diff.x, center.x),
//        fmaf(static_cast<float>(coord.y), diff.y, center.y), fmaf(static_cast<float>(coord.z), diff.z, center.z));
//    PKDParticle p;
//    p.dim = coord.dim;
//    p.pos = pos;
//    return p;
//}

struct QStackEntry {
    float t0, t1;
    unsigned int nodeID;
    box3f bounds;
};

MM_OPTIX_INTERSECTION_KERNEL(comp_treelets_intersect)
() {
    const int treeletID = optixGetPrimitiveIndex();
    const auto& self = getProgramData<QTreeletsGeoData>();
    const auto treelet = self.treeletBufferPtr[treeletID];

    const int begin = treelet.begin;
    const int size = treelet.end - begin;
    {
        float t0, t1;
        {
            auto const ray =
                Ray(optixGetWorldRayOrigin(), optixGetWorldRayDirection(), optixGetRayTmin(), optixGetRayTmax());
            if (!clipToBounds(ray, treelet.bounds, t0, t1))
                return;
        }
        auto const ray = Ray(optixGetWorldRayOrigin(), optixGetWorldRayDirection(), fmaxf(optixGetRayTmin(), t0),
            fminf(optixGetRayTmax(), t1));

        auto bounds = treelet.bounds;
        glm::vec3 _center;
        glm::vec3 _span;

        int nodeID = 0;
        float tmp_hit_t = ray.tmax;
        int tmp_hit_primID = -1;

        enum { STACK_DEPTH = 12 };
        QStackEntry stackBase[STACK_DEPTH];
        QStackEntry* stackPtr = stackBase;

        const int dir_sign[3] = {ray.direction.x < 0.f, ray.direction.y < 0.f, ray.direction.z < 0.f};
        const float org[3] = {ray.origin.x, ray.origin.y, ray.origin.z};
        const float rdir[3] = {
            (fabsf(ray.direction.x) <= 1e-8f) ? 1e8f : 1.f / ray.direction.x,
            (fabsf(ray.direction.y) <= 1e-8f) ? 1e8f : 1.f / ray.direction.y,
            (fabsf(ray.direction.z) <= 1e-8f) ? 1e8f : 1.f / ray.direction.z,
        };

        glm::vec3 tmp_hit_pos;

        while (1) {
            // while we have anything to traverse ...

            while (1) {
                // while we can go down
                _center = bounds.center();
                _span = bounds.span();

                const int particleID = nodeID + begin;
                const PKDParticle& particle = decode_coord(self.particleBufferPtr[particleID], _center, _span);
                int const dim = particle.dim;

                const float t_slab_lo = (particle.pos[dim] - self.radius - org[dim]) * rdir[dim];
                const float t_slab_hi = (particle.pos[dim] + self.radius - org[dim]) * rdir[dim];

                const float t_slab_nr = fminf(t_slab_lo, t_slab_hi);
                const float t_slab_fr = fmaxf(t_slab_lo, t_slab_hi);

                // -------------------------------------------------------
                // compute potential sphere interval, and intersect if necessary
                // -------------------------------------------------------
                const float sphere_t0 = fmaxf(t0, t_slab_nr);
                const float sphere_t1 = fminf(fminf(t_slab_fr, t1), tmp_hit_t);

                if (sphere_t0 < sphere_t1) {
                    if (intersectSphere(particle, self.radius, ray, tmp_hit_t)) {
                        tmp_hit_primID = particleID;
                        tmp_hit_pos = particle.pos;
                    }
                }

                // -------------------------------------------------------
                // compute near and far side intervals
                // -------------------------------------------------------
                const float nearSide_t0 = t0;
                const float nearSide_t1 = fminf(fminf(t_slab_fr, t1), tmp_hit_t);

                const float farSide_t0 = fmaxf(t0, t_slab_nr);
                const float farSide_t1 = fminf(t1, tmp_hit_t);

                // -------------------------------------------------------
                // logic
                // -------------------------------------------------------
                const int nearSide_nodeID = 2 * nodeID + 1 + dir_sign[dim];
                const int farSide_nodeID = 2 * nodeID + 2 - dir_sign[dim];

                const bool nearSide_valid = nearSide_nodeID < size;
                const bool farSide_valid = farSide_nodeID < size;

                const bool need_nearSide = nearSide_valid && nearSide_t0 < nearSide_t1;
                const bool need_farSide = farSide_valid && farSide_t0 < farSide_t1;

                // we have lB and rB
                // in case of dirSign == 1: near -> rB; far -> lB
                // in case of dirSign == 0: near -> lB; far -> rB

                if (!(need_nearSide || need_farSide))
                    break; // pop ...

                if (need_nearSide && need_farSide) {
                    stackPtr->nodeID = farSide_nodeID;
                    stackPtr->t0 = farSide_t0;
                    stackPtr->t1 = farSide_t1;
                    stackPtr->bounds = bounds;
                    if (dir_sign[dim]) {
                        // left
                        bounds.upper[dim] = particle.pos[dim] + self.radius;
                    } else {
                        // right
                        bounds.lower[dim] = particle.pos[dim] - self.radius;
                    }
                    ++stackPtr;

                    nodeID = nearSide_nodeID;
                    t0 = nearSide_t0;
                    t1 = nearSide_t1;
                    if (dir_sign[dim]) {
                        // right
                        bounds.lower[dim] = particle.pos[dim] - self.radius;
                    } else {
                        // left
                        bounds.upper[dim] = particle.pos[dim] + self.radius;
                    }
                    continue;
                }

                nodeID = need_nearSide ? nearSide_nodeID : farSide_nodeID;
                t0 = need_nearSide ? nearSide_t0 : farSide_t0;
                t1 = need_nearSide ? nearSide_t1 : farSide_t1;
                if (need_nearSide) {
                    if (dir_sign[dim]) {
                        // right
                        bounds.lower[dim] = particle.pos[dim] - self.radius;
                    } else {
                        // left
                        bounds.upper[dim] = particle.pos[dim] + self.radius;
                    }
                } else {
                    if (dir_sign[dim]) {
                        // left
                        bounds.upper[dim] = particle.pos[dim] + self.radius;
                    } else {
                        // right
                        bounds.lower[dim] = particle.pos[dim] - self.radius;
                    }
                }
            }
            // -------------------------------------------------------
            // pop
            // -------------------------------------------------------
            while (1) {
                if (stackPtr == stackBase) {
                    // can't pop any more - done.
                    if (tmp_hit_primID >= 0 && tmp_hit_t < ray.tmax) {
                        optixReportIntersection(tmp_hit_t, 0, tmp_hit_primID, __float_as_uint(tmp_hit_pos.x),
                            __float_as_uint(tmp_hit_pos.y), __float_as_uint(tmp_hit_pos.z));
                        //PerRayData& prd = getPerRayData<PerRayData>();
                        //prd.pos = tmp_hit_pos;
                    }
                    return;
                }
                --stackPtr;
                t0 = stackPtr->t0;
                t1 = stackPtr->t1;
                nodeID = stackPtr->nodeID;
                t1 = fminf(t1, tmp_hit_t);
                bounds = stackPtr->bounds;
                if (t1 <= t0)
                    continue;
                break;
            }
        }
    }
}

MM_OPTIX_CLOSESTHIT_KERNEL(comp_treelets_closesthit)
() {
    const unsigned int primID = optixGetAttribute_0();
    PerRayData& prd = getPerRayData<PerRayData>();

    const auto& self = getProgramData<QTreeletsGeoData>();

    prd.particleID = primID;
    //const PKDParticle& particle = self.particleBufferPtr[primID];
    //prd.pos = particle.pos;
    prd.pos.x = __uint_as_float(optixGetAttribute_1());
    prd.pos.y = __uint_as_float(optixGetAttribute_2());
    prd.pos.z = __uint_as_float(optixGetAttribute_3());
    glm::vec3 geo_col = glm::vec3(self.globalColor);
    if (self.hasColorData) {
        geo_col = glm::vec3(self.colorBufferPtr[primID]);
    }
    prd.albedo = geo_col;
    prd.t = optixGetRayTmax();
    set_depth(prd, optixGetRayTmax());
}


MM_OPTIX_CLOSESTHIT_KERNEL(comp_treelets_closesthit_occlusion)
() {
    optixSetPayload_0(1);
}

} // namespace device
} // namespace optix_hpg
} // namespace megamol
