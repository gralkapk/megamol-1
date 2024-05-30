#pragma once

#include <mutex>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "particle.h"
#include "pkd_utils.h"

namespace megamol::optix_hpg {
inline int arg_max2(glm::vec3 const& v) {
    int biggestDim = 0;
    for (int i = 1; i < 3; ++i)
        if ((v[i]) > (v[biggestDim]))
            biggestDim = i;
    return biggestDim;
}

inline size_t spkd_sort_partition(
    std::vector<device::PKDParticle>& particles, size_t begin, size_t end, device::box3f const& bounds, int& splitDim) {
    // -------------------------------------------------------
    // determine split pos
    // -------------------------------------------------------
    splitDim = arg_max2(bounds.upper - bounds.lower);

    float splitPos = (0.5f * (bounds.upper + bounds.lower))[splitDim];

    // -------------------------------------------------------
    // now partition ...
    // -------------------------------------------------------
    size_t mid = begin;
    size_t l = begin, r = (end - 1);
    // quicksort partition:
    while (l <= r) {
        while (l < r && particles[l].pos[splitDim] < splitPos)
            ++l;
        while (l < r && particles[r].pos[splitDim] >= splitPos)
            --r;
        if (l == r) {
            mid = l;
            break;
        }

        std::swap(particles[l], particles[r]);
    }

    // catch-all for extreme cases where all particles are on the same
    // spot, and can't be split:
    if (mid == begin || mid == end)
        mid = (begin + end) / 2;

    return mid;
}

template<typename MakeLeafLambda>
void partition_recurse(
    std::vector<device::PKDParticle>& particles, size_t begin, size_t end, const MakeLeafLambda& makeLeaf) {
    if (makeLeaf(begin, end, false))
        // could make into a leaf, done.
        return;

    // -------------------------------------------------------
    // parallel bounding box computation
    // -------------------------------------------------------
    device::box3f bounds;

    for (size_t idx = begin; idx < end; ++idx) {
        bounds.extend(particles[idx].pos);
    }

    int splitDim;
    auto mid = spkd_sort_partition(particles, begin, end, bounds, splitDim);

    // -------------------------------------------------------
    // and recurse ...
    // -------------------------------------------------------
    tbb::parallel_for(0, 2, [&](int side) {
        if (side)
            partition_recurse(particles, begin, mid, makeLeaf);
        else
            partition_recurse(particles, mid, end, makeLeaf);
    });
}

inline device::box3f extendBounds2(
    std::vector<device::PKDParticle> const& particles, size_t begin, size_t end, float radius) {
    device::box3f bounds;
    for (int64_t p_idx = begin; p_idx < end; ++p_idx) {
        auto const new_lower = particles[p_idx].pos - radius;
        auto const new_upper = particles[p_idx].pos + radius;
        bounds.extend(new_lower);
        bounds.extend(new_upper);
    }

    return bounds;
}

inline std::tuple<bool, std::unordered_set<unsigned char>, std::unordered_set<unsigned char>,
    std::unordered_set<unsigned char>>
prefix_consistency(
    std::vector<device::PKDParticle> const& particles, glm::vec3 const& lower, size_t begin, size_t end) {
    std::unordered_set<unsigned char> sx, sy, sz;
    byte_cast bc;
    bc.ui = 0;
    for (size_t i = begin; i < end; ++i) {
        auto const qpkd = encode_coord(particles[i].pos - lower, glm::vec3(), glm::vec3());
        bc.ui = qpkd.x;
        sx.insert(bc.parts.b);
        bc.ui = qpkd.y;
        sy.insert(bc.parts.b);
        bc.ui = qpkd.z;
        sz.insert(bc.parts.b);
    }
    return std::make_tuple(sx.size() <= device::spkd_array_size && sy.size() <= device::spkd_array_size &&
                               sz.size() <= device::spkd_array_size,
        sx, sy, sz);
}

inline std::vector<device::SPKDlet> partition_data(
    std::vector<device::PKDParticle>& particles, size_t maxSize, float radius) {
    std::mutex resultMutex;
    std::vector<device::SPKDlet> result;

    partition_recurse(particles, 0ULL, particles.size(), [&](size_t begin, size_t end, bool force) {
        /*bool makeLeaf() :*/
        const size_t size = end - begin;
        if (size > maxSize && !force)
            return false;

        auto [con, sx, sy, sz] = prefix_consistency(particles, glm::vec3(), begin, end);
        if (!con)
            return false;

        device::SPKDlet treelet;
        treelet.begin = begin;
        treelet.end = end;
        treelet.bounds = extendBounds2(particles, begin, end, radius);
        std::copy(sx.begin(), sx.end(), treelet.sx);
        std::copy(sy.begin(), sy.end(), treelet.sy);
        std::copy(sz.begin(), sz.end(), treelet.sz);

        std::lock_guard<std::mutex> lock(resultMutex);
        result.push_back(treelet);
        return true;
    });

    return std::move(result);
}

inline std::vector<device::SPKDlet> partition_data(std::vector<device::PKDParticle>& particles, size_t tbegin,
    size_t tend, glm::vec3 const& lower, size_t maxSize, float radius) {
    std::mutex resultMutex;
    std::vector<device::SPKDlet> result;

    partition_recurse(particles, tbegin, tend, [&](size_t begin, size_t end, bool force) {
        /*bool makeLeaf() :*/
        const size_t size = end - begin;
        if (size > maxSize && !force)
            return false;

        auto [con, sx, sy, sz] = prefix_consistency(particles, lower, begin, end);
        if (!con)
            return false;

        device::SPKDlet treelet;
        treelet.begin = begin;
        treelet.end = end;
        treelet.bounds = extendBounds2(particles, begin, end, radius);
        std::copy(sx.begin(), sx.end(), treelet.sx);
        std::copy(sy.begin(), sy.end(), treelet.sy);
        std::copy(sz.begin(), sz.end(), treelet.sz);

        std::lock_guard<std::mutex> lock(resultMutex);
        result.push_back(treelet);
        return true;
    });

    return std::move(result);
}


std::tuple<std::vector<glm::vec3>, std::vector<glm::vec3>, std::vector<glm::vec3>> compute_diffs(
    std::vector<device::SPKDlet> const& treelets, std::vector<device::SPKDParticle> const& sparticles,
    std::vector<device::PKDParticle> const& org_data, size_t gbegin, size_t gend) {
    std::vector<glm::vec3> diffs(gend - gbegin);
    std::vector<glm::vec3> ops(gend - gbegin);
    std::vector<glm::vec3> sps(gend - gbegin);
    tbb::parallel_for((size_t) 0, treelets.size(), [&](auto const tID) {
        //for (auto const& treelet : treelets) {
        auto const& treelet = treelets[tID];
        for (size_t i = treelet.begin; i < treelet.end; ++i) {
            glm::dvec3 const dpos = decode_spart(sparticles[i], treelet);
            glm::dvec3 const org_pos = org_data[i].pos;
            glm::dvec3 const diff = dpos - org_pos;
            diffs[i - gbegin] = diff;
            ops[i - gbegin] = org_pos;
            sps[i - gbegin] = dpos;
        }
        //}
    });
    return std::make_tuple(diffs, ops, sps);
}

} // namespace megamol::optix_hpg
