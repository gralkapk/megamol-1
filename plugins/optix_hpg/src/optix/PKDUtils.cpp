#include "PKDUtils.h"

#include <glm/glm.hpp>
#include <tbb/parallel_for.h>

namespace megamol::optix_hpg {
// BEGIN PKD

int arg_max(glm::vec3 const& v) {
    int biggestDim = 0;
    for (int i = 1; i < 3; ++i)
        if ((v[i]) > (v[biggestDim]))
            biggestDim = i;
    return biggestDim;
}

void recBuild(size_t /* root node */ P, device::PKDParticle* particle, size_t N, device::box3f bounds) {
    if (P >= N)
        return;

    int dim = arg_max(bounds.upper - bounds.lower);

    const size_t L = lChild(P);
    const size_t R = rChild(P);
    const bool lValid = (L < N);
    const bool rValid = (R < N);
    makeHeap(std::greater<float>(), L, particle, N, dim);
    makeHeap(std::less<float>(), R, particle, N, dim);

    if (rValid) {
        while (particle[L].pos[dim] > particle[R].pos[dim]) {
            std::swap(particle[L], particle[R]);
            trickle(std::greater<float>(), L, particle, N, dim);
            trickle(std::less<float>(), R, particle, N, dim);
        }
        if (particle[L].pos[dim] > particle[P].pos[dim]) {
            std::swap(particle[L], particle[P]);
            particle[L].dim = dim;
        } else if (particle[R].pos[dim] < particle[P].pos[dim]) {
            std::swap(particle[R], particle[P]);
            particle[R].dim = dim;
        } else
            /* nothing, root fits */;
    } else if (lValid) {
        if (particle[L].pos[dim] > particle[P].pos[dim]) {
            std::swap(particle[L], particle[P]);
            particle[L].dim = dim;
        }
    }

    device::box3f lBounds = bounds;
    device::box3f rBounds = bounds;
    lBounds.upper[dim] = rBounds.lower[dim] = particle[P].pos[dim];
    particle[P].dim = dim;

    tbb::parallel_for(0, 2, [&](int childID) {
        if (childID) {
            recBuild(L, particle, N, lBounds);
        } else {
            recBuild(R, particle, N, rBounds);
        }
    });
}

void makePKD(std::vector<device::PKDParticle>& particles, device::box3f bounds) {
    recBuild(/*node:*/ 0, particles.data(), particles.size(), bounds);
}

void makePKD(std::vector<device::PKDParticle>& particles, size_t begin, size_t end, device::box3f bounds) {
    recBuild(/*node:*/ 0, particles.data() + begin, end - begin, bounds);
}

// END PKD


// BEGIN TREELETS

size_t sort_partition(
    std::vector<device::PKDParticle>& particles, size_t begin, size_t end, device::box3f bounds, int& splitDim) {
    // -------------------------------------------------------
    // determine split pos
    // -------------------------------------------------------
    splitDim = arg_max(bounds.upper - bounds.lower);
    //float splitPos = bounds.center()[splitDim];
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

device::box3f extendBounds(std::vector<device::PKDParticle> const& particles, size_t begin, size_t end, float radius) {
    device::box3f bounds;
    for (int64_t p_idx = begin; p_idx < end; ++p_idx) {
        auto const new_lower = particles[p_idx].pos - radius;
        auto const new_upper = particles[p_idx].pos + radius;
        bounds.extend(new_lower);
        bounds.extend(new_upper);
    }

    return bounds;
}

std::vector<device::PKDlet> prePartition_inPlace(
    std::vector<device::PKDParticle>& particles, size_t maxSize, float radius) {
    std::mutex resultMutex;
    std::vector<device::PKDlet> result;

    partitionRecursively(particles, 0ULL, particles.size(), [&](size_t begin, size_t end, bool force) {
        /*bool makeLeaf() :*/
        const size_t size = end - begin;
        if (size > maxSize && !force)
            return false;

        device::PKDlet treelet;
        treelet.begin = begin;
        treelet.end = end;
        treelet.bounds = extendBounds(particles, begin, end, radius);

        std::lock_guard<std::mutex> lock(resultMutex);
        result.push_back(treelet);
        return true;
    });

    return std::move(result);
}

// END TREELETS

// BEGIN COMPRESSION

glm::uvec3 encode_coord(glm::vec3 const& pos, glm::vec3 const& center, glm::vec3 const& span) {
    constexpr unsigned int digits = 1023u;
    auto const dir = pos - center;
    auto const diff = span / static_cast<float>(digits);
    auto const coord = glm::uvec3(dir.x / diff.x, dir.y / diff.y, dir.z / diff.z);

    return coord;
}

glm::vec3 decode_coord(glm::uvec3 const& coord, glm::vec3 const& center, glm::vec3 const& span) {
    constexpr unsigned int digits = 1023u;
    auto const diff = span / static_cast<float>(digits);
    auto pos = glm::vec3(static_cast<float>(coord.x) * diff.x, static_cast<float>(coord.y) * diff.y,
        static_cast<float>(coord.z) * diff.z);
    pos = pos + center;

    return pos;
}

void convert(size_t P, device::PKDParticle* in_particle, device::QPKDParticle* out_particle, size_t N,
    device::box3f bounds,
    float radius, device::PKDParticle* out_decode) {
    if (P >= N)
        return;

    int const dim = in_particle[P].dim;

    auto const center = bounds.lower;
    auto const span = bounds.span();

    auto const coord = encode_coord(in_particle[P].pos, center, span);
    auto const pos = decode_coord(coord, center, span);
    if (out_decode) {
        out_decode[P].pos = pos;
    }

    out_particle[P].dim = dim;
    out_particle[P].x = coord.x;
    out_particle[P].y = coord.y;
    out_particle[P].z = coord.z;

    auto lBounds = bounds;
    auto rBounds = bounds;

    lBounds.upper[dim] = pos[dim] + radius;
    rBounds.lower[dim] = pos[dim] - radius;

    auto const L = lChild(P);
    auto const R = rChild(P);
    //const bool lValid = (L < N);
    //const bool rValid = (R < N);

    // TODO parallel
    convert(L, in_particle, out_particle, N, lBounds, radius, out_decode);
    convert(R, in_particle, out_particle, N, rBounds, radius, out_decode);
}

// END COMPRESSION
} // namespace megamol::optix_hpg
