#pragma once

#include <vector>
#include <tuple>
#include <mutex>

#include <glm/glm.hpp>

#include <geometry_calls/SimpleSphericalParticles.h>

namespace megamol::datatools {
inline int arg_max(glm::vec3 const& v) {
    int biggestDim = 0;
    for (int i = 1; i < 3; ++i)
        if ((v[i]) > (v[biggestDim]))
            biggestDim = i;
    return biggestDim;
}

typedef struct box3f {
    box3f()
            : lower{std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                  std::numeric_limits<float>::max()}
            , upper{std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(),
                  std::numeric_limits<float>::lowest()} {}
    ~box3f() = default;
    void extend(glm::vec3 const& pos) {
        lower = glm::vec3(fminf(pos.x, lower.x), fminf(pos.y, lower.y), fminf(pos.z, lower.z));
        upper = glm::vec3(fmaxf(pos.x, upper.x), fmaxf(pos.y, upper.y), fmaxf(pos.z, upper.z));
    }
    void extend(box3f const& box) {
        extend(box.lower);
        extend(box.upper);
    }
    glm::vec3 center() const {
        return (lower + upper) * 0.5f;
    }
    glm::vec3 span() const {
        return upper - lower;
    }
    glm::vec3 lower;
    glm::vec3 upper;
} box3f;

struct pkdlet {
    unsigned int begin, end;
    box3f bounds;
};

inline size_t lChild(size_t P) {
    return 2 * P + 1;
}
inline size_t rChild(size_t P) {
    return 2 * P + 2;
}

template<class Comp>
inline void trickle(
    const Comp& worse, size_t P, glm::vec3* particle, size_t N, int dim, glm::u8vec4* pack = nullptr) {
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
        if (rValid && worse(particle[R][dim], particle[L][dim]))
            C = R;

        if (!worse(particle[C][dim], particle[P][dim]))
            return;

        std::swap(particle[C], particle[P]);
        if (pack) {
            std::swap(pack[C], pack[P]);
        }
        P = C;
    }
}

template<class Comp>
inline void makeHeap(
    const Comp& comp, size_t P, glm::vec3* particle, size_t N, int dim, glm::u8vec4* pack = nullptr) {
    if (P >= N)
        return;
    const size_t L = lChild(P);
    const size_t R = rChild(P);
    makeHeap(comp, L, particle, N, dim, pack);
    makeHeap(comp, R, particle, N, dim, pack);
    trickle(comp, P, particle, N, dim, pack);
}

void recBuild(size_t P, glm::vec3* particle, size_t N, box3f bounds, glm::u8vec4* color = nullptr);

std::tuple<std::vector<glm::vec3>, std::vector<glm::u8vec4>, box3f> makePKD(
    geocalls::SimpleSphericalParticles& particles);

std::tuple<std::vector<glm::vec3>, std::vector<glm::u8vec4>, box3f> collectData(
    geocalls::SimpleSphericalParticles const& particles);

void makePKD(std::vector<glm::vec3>& position, std::vector<glm::u8vec4>& color, box3f const& bounds);

void makePKD(std::vector<glm::vec3>& position, std::vector<glm::u8vec4>& color, box3f const& bounds,
    unsigned int const begin, unsigned int const end);




inline size_t sort_partition(std::vector<glm::vec3>& particles, size_t begin, size_t end, box3f const& bounds, std::vector<glm::u8vec4>& color) {
    // -------------------------------------------------------
    // determine split pos
    // -------------------------------------------------------
    auto const span = bounds.span();
    auto const splitDim = arg_max(span);

    float splitPos = bounds.center()[splitDim];

    //float splitPos = (0.5f * (bounds.upper + bounds.lower))[splitDim];

    bool hasColor = !color.empty();

    // -------------------------------------------------------
    // now partition ...
    // -------------------------------------------------------
    size_t mid = begin;
    size_t l = begin, r = (end - 1);
    // quicksort partition:
    while (l <= r) {
        while (l < r && particles[l][splitDim] < splitPos)
            ++l;
        while (l < r && particles[r][splitDim] >= splitPos)
            --r;
        if (l == r) {
            mid = l;
            break;
        }

        std::swap(particles[l], particles[r]);
        if (hasColor) {
            std::swap(color[l], color[r]);
        }
    }

    // catch-all for extreme cases where all particles are on the same
    // spot, and can't be split:
    if (mid == begin || mid == end)
        mid = (begin + end) / 2;

    return mid;
}

template<typename MakeLeafLambda>
void partitionRecursively(std::vector<glm::vec3>& particles, size_t begin, size_t end, const MakeLeafLambda& makeLeaf, std::vector<glm::u8vec4>& color) {
    // -------------------------------------------------------
    // parallel bounding box computation
    // -------------------------------------------------------
    box3f bounds;

    for (size_t idx = begin; idx < end; ++idx) {
        bounds.extend(particles[idx]);
    }

    if (makeLeaf(begin, end, false, bounds))
        // could make into a leaf, done.
        return;

    auto mid = sort_partition(particles, begin, end, bounds, color);

    // -------------------------------------------------------
    // and recurse ...
    // -------------------------------------------------------
    tbb::parallel_for(0, 2, [&](int side) {
        if (side)
            partitionRecursively(particles, begin, mid, makeLeaf, color);
        else
            partitionRecursively(particles, mid, end, makeLeaf, color);
    });
}

std::vector<pkdlet> prePartition_inPlace(std::vector<glm::vec3>& particles, size_t maxSize, float radius, std::vector<glm::u8vec4>& color);
} // namespace megamol::datatools
