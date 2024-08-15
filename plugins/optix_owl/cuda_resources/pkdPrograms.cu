#include <cuda_runtime.h>
#include <optix_device.h>

#include <owl/common/math/box.h>
#include <owl/common/math/vec.h>
#include <owl/owl_device.h>

#include "perraydata.h"
#include "pkd.h"
#include "intersect.h"

namespace megamol {
namespace optix_owl {
namespace device {
using namespace owl::common;

struct StackEntry {
    float t0, t1;
    unsigned int nodeID;
};

OPTIX_INTERSECT_PROGRAM(pkd_intersect)() {
    const auto& self = owl::getProgramData<PKDGeomData>();

    owl::Ray ray(optixGetWorldRayOrigin(), optixGetWorldRayDirection(), optixGetRayTmin(), optixGetRayTmax());

    float t0, t1;
    if (!clipToBounds(ray, self.worldBounds, t0, t1))
        return;


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
    float const particleRadius = self.particleRadius;

    while (1) {
        // while we have anything to traverse ...

        while (1) {
            // while we can go down
            const Particle particle = self.particleBuffer[nodeID];
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
                if (intersectSphere(particle, particleRadius, ray, tmp_hit_t)) {
                    tmp_hit_primID = nodeID;
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

            const bool nearSide_valid = nearSide_nodeID < numParticles;
            const bool farSide_valid = farSide_nodeID < numParticles;

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
            t1 = fminf(stackPtr->t1, tmp_hit_t);
            nodeID = stackPtr->nodeID;
            if (t1 <= t0)
                continue;
            break;
        }
    }
}

OPTIX_CLOSEST_HIT_PROGRAM(pkd_ch)() {
    PerRayData& prd = owl::getPRD<PerRayData>();
    const auto& self = owl::getProgramData<PKDGeomData>();
    prd.particleID = optixGetAttribute_0();
    prd.t = optixGetRayTmax();
    prd.pos = self.particleBuffer[optixGetAttribute_0()].pos;
}

OPTIX_BOUNDS_PROGRAM(pkd_bounds)(const void* geomData, box3f& primBounds, const int primID) {
    auto& self = *(const PKDGeomData*) geomData;
    primBounds = self.worldBounds;
}
} // namespace device
} // namespace optix_owl
} // namespace megamol
