#pragma once

#include <vector>

#include "particle.h"
#include "utils_host.h"

namespace megamol::optix_hpg {
// BEGIN PKD
void makePKD(std::vector<device::PKDParticle>& particles, box3f bounds);

void makePKD(std::vector<device::PKDParticle>& particles, size_t begin, size_t end, box3f bounds);

std::vector<device::PKDlet> prePartition_inPlace(
    std::vector<device::PKDParticle>& particles, size_t maxSize, float radius);
// END PKD

// BEGIN TREELETS
size_t sort_partition(
    std::vector<device::PKDParticle>& particles, size_t begin, size_t end, box3f bounds, int& splitDim);

template<typename MakeLeafLambda>
void partitionRecursively(
    std::vector<device::PKDParticle>& particles, size_t begin, size_t end, const MakeLeafLambda& makeLeaf) {
    if (makeLeaf(begin, end, false))
        // could make into a leaf, done.
        return;

    // -------------------------------------------------------
    // parallel bounding box computation
    // -------------------------------------------------------
    box3f bounds;

    for (size_t idx = begin; idx < end; ++idx) {
        bounds.extend(particles[idx].pos);
    }

#if 0
    std::mutex boundsMutex;

    const size_t blockSize = 32 * 1024;
    const size_t numTasks = end - begin;
    const size_t numBlocks = (numTasks + blockSize - 1) / blockSize;
    /*parallel_for(numBlocks, [&](size_t blockID) {
        size_t block_begin = begin + blockID * blockSize;
        taskFunction(block_begin, std::min(block_begin + blockSize, end));
    });*/

    auto task = [&](size_t blockBegin, size_t blockEnd) {
        //box3f blockBounds;
        /*for (size_t i = blockBegin; i < blockEnd; i++)
            blockBounds.extend(model->particles[i].pos);*/
        auto const blockBounds = extendBounds(particles, blockBegin, blockEnd, 0.f);
        std::lock_guard<std::mutex> lock(boundsMutex);
        bounds.extend(blockBounds);
    };
    tbb::parallel_for((size_t) 0, numBlocks, [&](size_t blockID) {
        size_t block_begin = begin + blockID * blockSize;
        task(block_begin, std::min(block_begin + blockSize, end));
    });

    //parallel_for_blocked(begin, end, 32 * 1024, [&](size_t blockBegin, size_t blockEnd) {
    //    box3f blockBounds;
    //    /*for (size_t i = blockBegin; i < blockEnd; i++)
    //        blockBounds.extend(model->particles[i].pos);*/
    //    std::tie(blockBounds.lower, blockBounds.upper) = extendBounds(particles, blockBegin, blockEnd, 0.f);
    //    std::lock_guard<std::mutex> lock(boundsMutex);
    //    /*bounds.extend(blockBounds);*/
    //    bounds.lower.x = std::min(bounds.lower.x, blockBounds.lower.x);
    //    bounds.lower.y = std::min(bounds.lower.y, blockBounds.lower.y);
    //    bounds.lower.z = std::min(bounds.lower.z, blockBounds.lower.z);
    //    bounds.upper.x = std::max(bounds.upper.x, blockBounds.upper.x);
    //    bounds.upper.y = std::max(bounds.upper.y, blockBounds.upper.y);
    //    bounds.upper.z = std::max(bounds.upper.z, blockBounds.upper.z);
    //});
#endif

    int splitDim;
    auto mid = sort_partition(particles, begin, end, bounds, splitDim);

    // -------------------------------------------------------
    // and recurse ...
    // -------------------------------------------------------
    tbb::parallel_for(0, 2, [&](int side) {
        if (side)
            partitionRecursively(particles, begin, mid, makeLeaf);
        else
            partitionRecursively(particles, mid, end, makeLeaf);
    });
}
// END TREELETS

// BEGIN COMPRESS
void convert(size_t P, device::PKDParticle* in_particle, device::QPKDParticle* out_particle, size_t N, box3f bounds,
    float radius);
// END COMPRESS

inline size_t parent(size_t C) {
    if (C == 0)
        return 0;
    return (C - 1) / 2;
}

inline size_t lChild(size_t P) {
    return 2 * P + 1;
}
inline size_t rChild(size_t P) {
    return 2 * P + 2;
}

template<class Comp>
inline void trickle(const Comp& worse, size_t P, device::PKDParticle* particle, size_t N, int dim) {
    if (P >= N)
        return;

    while (1) {
        const size_t L = lChild(P);
        const size_t R = rChild(P);
        const bool lValid = (L < N);
        const bool rValid = (R < N);

        if (!lValid)
            return;
        size_t C = L;
        if (rValid && worse(particle[R].pos[dim], particle[L].pos[dim]))
            C = R;

        if (!worse(particle[C].pos[dim], particle[P].pos[dim]))
            return;

        std::swap(particle[C], particle[P]);
        P = C;
    }
}

template<class Comp>
inline void makeHeap(const Comp& comp, size_t P, device::PKDParticle* particle, size_t N, int dim) {
    if (P >= N)
        return;
    const size_t L = lChild(P);
    const size_t R = rChild(P);
    makeHeap(comp, L, particle, N, dim);
    makeHeap(comp, R, particle, N, dim);
    trickle(comp, P, particle, N, dim);
}
} // namespace megamol::optix_hpg
