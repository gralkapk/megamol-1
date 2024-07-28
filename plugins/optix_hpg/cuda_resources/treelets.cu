#include "pkd.h"

#include "optix/utils_device.h"
#include "pkd_utils.h"
#include "ptx_util.h"

#include "datatools/PKDUtils.h"

namespace megamol {
namespace optix_hpg {
namespace device {

MM_OPTIX_INTERSECTION_KERNEL(treelets_intersect_flat)() {
    const int treeletID = optixGetPrimitiveIndex();

    const auto& self = getProgramData<TreeletsGeoData>();

    const auto treelet = self.treeletBufferPtr[treeletID];

    auto const ray = Ray(optixGetWorldRayOrigin(), optixGetWorldRayDirection(), optixGetRayTmin(), optixGetRayTmax());
    float t0, t1;
    if (!clipToBounds(ray, treelet.bounds, t0, t1))
        return;

    float tmp_hit_t = ray.tmax;
    int tmp_hit_primID = -1;

    for (unsigned int i = treelet.begin; i < treelet.end; ++i) {
        const glm::vec3& particle = self.particleBufferPtr[i];
        if (intersectSphere(particle, self.radius, ray, tmp_hit_t)) {
            tmp_hit_primID = i;
        }
    }

    if (tmp_hit_primID >= 0 && tmp_hit_t < ray.tmax) {
        optixReportIntersection(tmp_hit_t, 0, tmp_hit_primID);
    }
}


MM_OPTIX_INTERSECTION_KERNEL(treelets_intersect)
() {
    const int treeletID = optixGetPrimitiveIndex();
    const auto& self = getProgramData<TreeletsGeoData>();
    const auto treelet = self.treeletBufferPtr[treeletID];

    const int begin = treelet.begin;
    const int size = treelet.end - begin;
    {
        float t0, t1;
        //{
        auto const ray =
            Ray(optixGetWorldRayOrigin(), optixGetWorldRayDirection(), optixGetRayTmin(), optixGetRayTmax());
        if (!clipToBounds(ray, treelet.bounds, t0, t1))
            return;
        //}
        /*auto const ray = Ray(optixGetWorldRayOrigin(), optixGetWorldRayDirection(), fmaxf(optixGetRayTmin(), t0),
            fminf(optixGetRayTmax(), t1));*/

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
                const glm::vec3& particle = self.particleBufferPtr[particleID];
                int const dim = getDim(particle);

                const float t_slab_lo = (particle[dim] - self.radius - org[dim]) * rdir[dim];
                const float t_slab_hi = (particle[dim] + self.radius - org[dim]) * rdir[dim];

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
    const glm::vec3& particle = self.particleBufferPtr[primID];
    prd.pos = particle;
    glm::vec3 geo_col = glm::vec3(self.globalColor) / 255.f;
    if (self.hasColorData) {
        geo_col = glm::vec3(self.colorBufferPtr[primID]) / 255.f;
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

} // namespace device
} // namespace optix_hpg
} // namespace megamol
