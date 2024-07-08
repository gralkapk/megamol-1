#pragma once

#include <functional>
#include <vector>

#include "particle.h"
#include "utils_host.h"

namespace megamol::optix_hpg {
// BEGIN PKD
template<typename PackType>
void recBuild(
    size_t /* root node */ P, device::PKDParticle* particle, size_t N, device::box3f bounds, PackType* pack = nullptr) {
    if (P >= N)
        return;

    int dim = arg_max(bounds.upper - bounds.lower);

    const size_t L = lChild(P);
    const size_t R = rChild(P);
    const bool lValid = (L < N);
    const bool rValid = (R < N);
    makeHeap(std::greater<float>(), L, particle, N, dim, pack);
    makeHeap(std::less<float>(), R, particle, N, dim, pack);

    if (rValid) {
        while (particle[L].pos[dim] > particle[R].pos[dim]) {
            std::swap(particle[L], particle[R]);
            if (pack) {
                std::swap(pack[L], pack[R]);
            }
            trickle(std::greater<float>(), L, particle, N, dim, pack);
            trickle(std::less<float>(), R, particle, N, dim, pack);
        }
        if (particle[L].pos[dim] > particle[P].pos[dim]) {
            std::swap(particle[L], particle[P]);
            if (pack) {
                std::swap(pack[L], pack[P]);
                pack[L].dim = dim;
            }
            particle[L].dim = dim;
        } else if (particle[R].pos[dim] < particle[P].pos[dim]) {
            std::swap(particle[R], particle[P]);
            if (pack) {
                std::swap(pack[R], pack[P]);
                pack[R].dim = dim;
            }
            particle[R].dim = dim;
        } else
            /* nothing, root fits */;
    } else if (lValid) {
        if (particle[L].pos[dim] > particle[P].pos[dim]) {
            std::swap(particle[L], particle[P]);
            if (pack) {
                std::swap(pack[L], pack[P]);
                pack[L].dim = dim;
            }
            particle[L].dim = dim;
        }
    }

    device::box3f lBounds = bounds;
    device::box3f rBounds = bounds;
    lBounds.upper[dim] = rBounds.lower[dim] = particle[P].pos[dim];
    particle[P].dim = dim;

    tbb::parallel_for(0, 2, [&](int childID) {
        if (childID) {
            recBuild(L, particle, N, lBounds, pack);
        } else {
            recBuild(R, particle, N, rBounds, pack);
        }
    });
}

void makePKD(std::vector<device::PKDParticle>& particles, device::box3f bounds);

template<typename PackType>
void makePKD(std::vector<device::PKDParticle>& particles, size_t begin, size_t end, device::box3f bounds,
    PackType* pack = nullptr) {
    recBuild(/*node:*/ 0, particles.data() + begin, end - begin, bounds, pack);
}

void makePKD(std::vector<device::SPKDParticle>& particles, device::SPKDlet const& treelet, size_t begin);

void makePKD(std::vector<device::C2PKDParticle>& particles, device::C2PKDlet const& treelet,
             device::box3f const& global_bounds, device::MortonConfig const& config);

std::vector<device::PKDlet> prePartition_inPlace(std::vector<device::PKDParticle>& particles, size_t maxSize,
    float radius, std::function<bool(device::box3f const&)> add_cond = nullptr);

std::vector<device::PKDlet> prePartition_inPlace(std::vector<device::PKDParticle>& particles, size_t begin, size_t end,
    size_t maxSize, float radius, std::function<bool(device::box3f const&)> add_cond = nullptr);
// END PKD

// BEGIN TREELETS
inline int arg_max(glm::vec3 const& v) {
    int biggestDim = 0;
    for (int i = 1; i < 3; ++i)
        if ((v[i]) > (v[biggestDim]))
            biggestDim = i;
    return biggestDim;
}

inline int arg_max(glm::uvec3 const& v) {
    int biggestDim = 0;
    for (int i = 1; i < 3; ++i)
        if ((v[i]) > (v[biggestDim]))
            biggestDim = i;
    return biggestDim;
}

template<typename PType, typename BType>
size_t sort_partition(std::vector<PType>& particles, size_t begin, size_t end, BType const& bounds, int& splitDim,
    bool median, device::PKDParticle* pack = nullptr) {
    // -------------------------------------------------------
    // determine split pos
    // -------------------------------------------------------
    auto const span = bounds.span();
    splitDim = arg_max(span);

    float splitPos = bounds.center()[splitDim];

    if (median) {
        std::nth_element(particles.begin() + begin, particles.begin() + begin + ((end - begin) / 2),
            particles.begin() + end,
            [&](auto const& lhs, auto const& rhs) { return lhs.pos[splitDim] < rhs.pos[splitDim]; });
        splitPos = (particles.begin() + begin + ((end - begin) / 2))->pos[splitDim];
    }

    //float splitPos = (0.5f * (bounds.upper + bounds.lower))[splitDim];

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
        if (pack) {
            std::swap(pack[l], pack[r]);
        }
    }

    // catch-all for extreme cases where all particles are on the same
    // spot, and can't be split:
    if (mid == begin || mid == end)
        mid = (begin + end) / 2;

    return mid;
}

template<typename PType, typename BType, typename MakeLeafLambda>
void partitionRecursively(std::vector<PType>& particles, size_t begin, size_t end, const MakeLeafLambda& makeLeaf,
    device::PKDParticle* pack = nullptr, bool median = false, bool skip = false) {
    // -------------------------------------------------------
    // parallel bounding box computation
    // -------------------------------------------------------
    BType bounds;

    for (size_t idx = begin; idx < end; ++idx) {
        bounds.extend(particles[idx].pos);
    }

    if (!skip) {
        if (makeLeaf(begin, end, false, bounds))
            // could make into a leaf, done.
            return;
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
    auto mid = sort_partition(particles, begin, end, bounds, splitDim, median, pack);

    // -------------------------------------------------------
    // and recurse ...
    // -------------------------------------------------------
    tbb::parallel_for(0, 2, [&](int side) {
        if (side)
            partitionRecursively<PType, BType>(particles, begin, mid, makeLeaf, pack, median);
        else
            partitionRecursively<PType, BType>(particles, mid, end, makeLeaf, pack, median);
    });
}
// END TREELETS

// BEGIN COMPRESS

std::vector<std::pair<size_t, size_t>> gridify(
    std::vector<device::PKDParticle>& data, glm::vec3 const& lower, glm::vec3 const& upper);
device::box3f extendBounds(std::vector<device::PKDParticle> const& particles, size_t begin, size_t end, float radius);
std::tuple<std::vector<device::PKDlet>, std::vector<std::pair<unsigned int, device::QPKDParticle>>> makeSpartition(
    std::vector<device::QPKDParticle> const& data, size_t begin, size_t end, float radius);
//std::tuple<std::vector<device::SPKDlet>, std::vector<device::SPKDParticle>> slice_qparticles(
//    std::vector<device::PKDlet> const& treelets,
//    std::vector<std::pair<unsigned int, device::QPKDParticle>> const& particles,
//    std::vector<device::PKDParticle> const& org_data, size_t begin, size_t end, float radius);
//std::vector<glm::vec3> compute_diffs(std::vector<device::SPKDlet> const& treelets,
//    std::vector<device::SPKDParticle> const& sparticles,
//    std::vector<std::pair<unsigned int, device::QPKDParticle>> const& qparticles,
//    std::vector<device::PKDParticle> const& org_data, size_t begin, size_t end, glm::vec3 const& lower);

void convert(size_t P, device::PKDParticle* in_particle, device::QPKDParticle* out_particle, size_t N,
    device::box3f bounds, float radius, device::PKDParticle* out_decode = nullptr, glm::uvec3* out_coord = nullptr);
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

template<class Comp, typename PackType>
inline void trickle(
    const Comp& worse, size_t P, device::PKDParticle* particle, size_t N, int dim, PackType* pack = nullptr) {
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
        if (pack) {
            std::swap(pack[C], pack[P]);
        }
        P = C;
    }
}

template<class Comp, typename PackType>
inline void makeHeap(
    const Comp& comp, size_t P, device::PKDParticle* particle, size_t N, int dim, PackType* pack = nullptr) {
    if (P >= N)
        return;
    const size_t L = lChild(P);
    const size_t R = rChild(P);
    makeHeap(comp, L, particle, N, dim, pack);
    makeHeap(comp, R, particle, N, dim, pack);
    trickle(comp, P, particle, N, dim, pack);
}


template<class Comp>
inline void trickle(
    const Comp& worse, size_t P, device::SPKDParticle* particle, size_t N, int dim, device::SPKDlet const& treelet) {
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
        if (rValid && worse(decode_spart(particle[R], treelet)[dim], decode_spart(particle[L], treelet)[dim]))
            C = R;

        if (!worse(decode_spart(particle[C], treelet)[dim], decode_spart(particle[P], treelet)[dim]))
            return;

        std::swap(particle[C], particle[P]);
        P = C;
    }
}

template<class Comp>
inline void makeHeap(
    const Comp& comp, size_t P, device::SPKDParticle* particle, size_t N, int dim, device::SPKDlet const& treelet) {
    if (P >= N)
        return;
    const size_t L = lChild(P);
    const size_t R = rChild(P);
    makeHeap(comp, L, particle, N, dim, treelet);
    makeHeap(comp, R, particle, N, dim, treelet);
    trickle(comp, P, particle, N, dim, treelet);
}


template<class Comp>
inline void trickle(
    const Comp& worse, size_t P, device::C2PKDParticle* particle, size_t N, int dim, device::C2PKDlet const& treelet, device::box3f const& global_bounds, device::MortonConfig const& config) {
    if (P >= N)
        return;

    auto const span = global_bounds.span();
    auto const lower = global_bounds.lower;

    while (1) {
        const size_t L = lChild(P);
        const size_t R = rChild(P);
        const bool lValid = (L < N);
        const bool rValid = (R < N);

        if (!lValid)
            return;
        size_t C = L;
        if (rValid &&
            worse((particle[R].from(
                      treelet.prefix, span, lower))[dim],
                          (particle[L].from(treelet.prefix, span, lower))[dim]))
            C = R;

        if (!worse((particle[C].from(
                       treelet.prefix, span, lower))[dim],
                (particle[P].from(
                    treelet.prefix, span, lower))[dim]))
            return;

        std::swap(particle[C], particle[P]);
        P = C;
    }
}

template<class Comp>
inline void makeHeap(
    const Comp& comp, size_t P, device::C2PKDParticle* particle, size_t N, int dim, device::C2PKDlet const& treelet, device::box3f const& global_bounds, device::MortonConfig const& config) {
    if (P >= N)
        return;
    const size_t L = lChild(P);
    const size_t R = rChild(P);
    makeHeap(comp, L, particle, N, dim, treelet, global_bounds, config);
    makeHeap(comp, R, particle, N, dim, treelet, global_bounds, config);
    trickle(comp, P, particle, N, dim, treelet, global_bounds, config);
}
} // namespace megamol::optix_hpg