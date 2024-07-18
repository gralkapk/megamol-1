#pragma once

#include <vector>
#include <tuple>

#include <glm/glm.hpp>

#include <geometry_calls/SimpleSphericalParticles.h>

namespace megamol::datatools {
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

std::tuple<std::vector<glm::vec3>, std::vector<glm::u8vec4>, box3f> makePKD(geocalls::SimpleSphericalParticles& particles);
} // namespace megamol::datatools
