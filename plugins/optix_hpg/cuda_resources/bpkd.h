#pragma once

struct BPKDStackEntry {
    float t0, t1;
    unsigned int nodeID;
    box3f refBox;
};


MM_OPTIX_INTERSECTION_KERNEL(treelet_intersect_bpkd)() {
    const int treeletID = optixGetPrimitiveIndex();
    const auto& self = getProgramData<BTreeletsGeoData>();
    const auto treelet = self.treeletBufferPtr[treeletID];

    auto const ray = Ray(optixGetWorldRayOrigin(), optixGetWorldRayDirection(), optixGetRayTmin(), optixGetRayTmax());

    const int begin = treelet.begin;
    const int size = treelet.end - begin;
    {
        float t0, t1;
        if (!clipToBounds(ray, treelet.bounds, t0, t1))
            return;


        int nodeID = 0;
        float tmp_hit_t = ray.tmax;
        int tmp_hit_primID = -1;

        enum { STACK_DEPTH = 12 };

        BPKDStackEntry stackBase[STACK_DEPTH];
        BPKDStackEntry* stackPtr = stackBase;

        glm::vec3 pos;
        int dim;

        box3f refBox = treelet.bounds;

        glm::vec3 tmp_hit_pos;

        while (1) {
            // while we have anything to traverse ...

            while (1) {
                // while we can go down

                const int particleID = nodeID + begin;

                {
                    auto const& bpart = self.particleBufferPtr[particleID];
                    dim = bpart.dim;
                    pos = bpart.from(refBox.span(), refBox.lower);
                }

                const float t_slab_lo =
                    (pos[dim] - self.radius - ray.origin[dim]) / ray.direction[dim] - t_compensate(refBox.span()[dim]);
                const float t_slab_hi =
                    (pos[dim] + self.radius - ray.origin[dim]) / ray.direction[dim] + t_compensate(refBox.span()[dim]);

                const float t_slab_nr = fminf(t_slab_lo, t_slab_hi);
                const float t_slab_fr = fmaxf(t_slab_lo, t_slab_hi);

                // -------------------------------------------------------
                // compute potential sphere interval, and intersect if necessary
                // -------------------------------------------------------
                /*const float sphere_t0 = fmaxf(t0, t_slab_nr);
                const float sphere_t1 = fminf(fminf(t_slab_fr, t1), tmp_hit_t);*/

                //if (sphere_t0 < sphere_t1) {
                if (intersectSphere(pos, self.radius, ray, tmp_hit_t)) {
                    tmp_hit_primID = particleID;

                    tmp_hit_pos = pos;
                }
                //}

                // -------------------------------------------------------
                // compute near and far side intervals
                // -------------------------------------------------------
                const float nearSide_t0 = t0;
                const float nearSide_t1 = fminf(fminf(t_slab_fr, t1), tmp_hit_t - t_compensate(refBox.span()[dim]));

                const float farSide_t0 = fmaxf(t0, t_slab_nr);
                const float farSide_t1 = fminf(t1, tmp_hit_t + t_compensate(refBox.span()[dim]));

                // -------------------------------------------------------
                // logic
                // -------------------------------------------------------
                const int nearSide_nodeID = 2 * nodeID + 1 + (ray.direction[dim] < 0.f);
                const int farSide_nodeID = 2 * nodeID + 2 - (ray.direction[dim] < 0.f);

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
                    stackPtr->refBox = rightBounds(refBox, pos[dim], self.radius, dim);

                    ++stackPtr;

                    nodeID = nearSide_nodeID;
                    t0 = nearSide_t0;
                    t1 = nearSide_t1;
                    refBox = leftBounds(refBox, pos[dim], self.radius, dim);

                    continue;
                }

                nodeID = need_nearSide ? nearSide_nodeID : farSide_nodeID;
                t0 = need_nearSide ? nearSide_t0 : farSide_t0;
                t1 = need_nearSide ? nearSide_t1 : farSide_t1;
                refBox = need_nearSide ? leftBounds(refBox, pos[dim], self.radius, dim)
                                       : rightBounds(refBox, pos[dim], self.radius, dim);
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
                    }
                    return;
                }
                --stackPtr;
                //getEntry(stackPtr, nodeID, t0, t1);
                t0 = stackPtr->t0;
                t1 = stackPtr->t1;
                nodeID = stackPtr->nodeID;
                t1 = min(t1, tmp_hit_t);

                refBox = stackPtr->refBox;

                if (t1 <= t0)
                    continue;
                break;
            }
        }
    }
}


MM_OPTIX_CLOSESTHIT_KERNEL(bpkd_treelets_closesthit)
() {
    const unsigned int primID = optixGetAttribute_0();
    PerRayData& prd = getPerRayData<PerRayData>();

    const auto& self = getProgramData<BTreeletsGeoData>();

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


MM_OPTIX_CLOSESTHIT_KERNEL(bpkd_treelets_closesthit_occlusion)
() {
    optixSetPayload_0(1);
}
