#include "datatools/PKD.h"

#include <type_traits>

#include <omp.h>
#include <tbb/parallel_for.h>

namespace megamol::datatools {
void set_dim(glm::vec3& p, int dim) {
    int& pxAsInt = reinterpret_cast<int&>(p.x);
    pxAsInt = (pxAsInt & ~3) | dim;
}

void recBuild(size_t /* root node */ P, glm::vec3* particle, size_t N, box3f bounds, glm::u8vec4* pack) {
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
        while (particle[L][dim] > particle[R][dim]) {
            std::swap(particle[L], particle[R]);
            if (pack) {
                std::swap(pack[L], pack[R]);
            }
            trickle(std::greater<float>(), L, particle, N, dim, pack);
            trickle(std::less<float>(), R, particle, N, dim, pack);
        }
        if (particle[L][dim] > particle[P][dim]) {
            std::swap(particle[L], particle[P]);
            if (pack) {
                std::swap(pack[L], pack[P]);
            }
            set_dim(particle[L], dim);
        } else if (particle[R][dim] < particle[P][dim]) {
            std::swap(particle[R], particle[P]);
            if (pack) {
                std::swap(pack[R], pack[P]);
            }
            set_dim(particle[R], dim);
        } else
            /* nothing, root fits */;
    } else if (lValid) {
        if (particle[L][dim] > particle[P][dim]) {
            std::swap(particle[L], particle[P]);
            if (pack) {
                std::swap(pack[L], pack[P]);
            }
            set_dim(particle[L], dim);
        }
    }

    box3f lBounds = bounds;
    box3f rBounds = bounds;
    lBounds.upper[dim] = rBounds.lower[dim] = particle[P][dim];
    set_dim(particle[P], dim);

    tbb::parallel_for(0, 2, [&](int childID) {
        if (childID) {
            recBuild(L, particle, N, lBounds, pack);
        } else {
            recBuild(R, particle, N, rBounds, pack);
        }
    });
}

bool has_global_color(geocalls::SimpleSphericalParticles::ColourDataType const& type) {
    return type == geocalls::SimpleSphericalParticles::ColourDataType::COLDATA_DOUBLE_I ||
           type == geocalls::SimpleSphericalParticles::ColourDataType::COLDATA_FLOAT_I ||
           type == geocalls::SimpleSphericalParticles::ColourDataType::COLDATA_NONE;
}

std::tuple<std::vector<glm::vec3>, std::vector<glm::u8vec4>, box3f> makePKD(geocalls::SimpleSphericalParticles& particles) {
    auto const p_count = particles.GetCount();
    std::vector<glm::vec3> position(p_count);
    
    auto const x_acc = particles.GetParticleStore().GetXAcc();
    auto const y_acc = particles.GetParticleStore().GetYAcc();
    auto const z_acc = particles.GetParticleStore().GetZAcc();

    auto const t_num = omp_get_max_threads();
    auto const t_total_num = std::min(static_cast<uint64_t>(t_num), p_count);
    std::vector<box3f> boxes(t_total_num);

#pragma omp parallel for num_threads(t_total_num)
    for (int64_t i = 0; i < p_count; ++i) {
        position[i].x = x_acc->Get_f(i);
        position[i].y = y_acc->Get_f(i);
        position[i].z = z_acc->Get_f(i);

        boxes[omp_get_thread_num()].extend(position[i]);
    }

    box3f bounds;
    for (int i = 0; i < t_total_num; ++i) {
        bounds.extend(boxes[i]);
    }

    if (!has_global_color(particles.GetColourDataType())) {
        std::vector<glm::u8vec4> color(p_count);

        auto const cr_acc = particles.GetParticleStore().GetCRAcc();
        auto const cg_acc = particles.GetParticleStore().GetCGAcc();
        auto const cb_acc = particles.GetParticleStore().GetCBAcc();
        auto const ca_acc = particles.GetParticleStore().GetCAAcc();

#pragma omp parallel for num_threads(t_total_num)
        for (int64_t i = 0; i < p_count; ++i) {
            color[i].r = cr_acc->Get_u8(i);
            color[i].g = cg_acc->Get_u8(i);
            color[i].b = cb_acc->Get_u8(i);
            color[i].a = ca_acc->Get_u8(i);
        }

        recBuild(0, position.data(), position.size(), bounds, color.data());
        return std::make_tuple(position, color, bounds);
    } else {
        recBuild(0, position.data(), position.size(), bounds);
        return std::make_tuple(position, std::vector<glm::u8vec4>(), bounds);
    }
}

std::tuple<std::vector<glm::vec3>, std::vector<glm::u8vec4>, box3f> collectData(
    geocalls::SimpleSphericalParticles const& particles) {
    auto const p_count = particles.GetCount();
    std::vector<glm::vec3> position(p_count);

    auto const x_acc = particles.GetParticleStore().GetXAcc();
    auto const y_acc = particles.GetParticleStore().GetYAcc();
    auto const z_acc = particles.GetParticleStore().GetZAcc();

    auto const t_num = omp_get_max_threads();
    auto const t_total_num = std::min(static_cast<uint64_t>(t_num), p_count);
    std::vector<box3f> boxes(t_total_num);

#pragma omp parallel for num_threads(t_total_num)
    for (int64_t i = 0; i < p_count; ++i) {
        position[i].x = x_acc->Get_f(i);
        position[i].y = y_acc->Get_f(i);
        position[i].z = z_acc->Get_f(i);

        boxes[omp_get_thread_num()].extend(position[i]);
    }

    box3f bounds;
    for (int i = 0; i < t_total_num; ++i) {
        bounds.extend(boxes[i]);
    }

    if (!has_global_color(particles.GetColourDataType())) {
        std::vector<glm::u8vec4> color(p_count);

        auto const cr_acc = particles.GetParticleStore().GetCRAcc();
        auto const cg_acc = particles.GetParticleStore().GetCGAcc();
        auto const cb_acc = particles.GetParticleStore().GetCBAcc();
        auto const ca_acc = particles.GetParticleStore().GetCAAcc();

#pragma omp parallel for num_threads(t_total_num)
        for (int64_t i = 0; i < p_count; ++i) {
            color[i].r = cr_acc->Get_u8(i);
            color[i].g = cg_acc->Get_u8(i);
            color[i].b = cb_acc->Get_u8(i);
            color[i].a = ca_acc->Get_u8(i);
        }
        return std::make_tuple(position, color, bounds);
    } else {
        return std::make_tuple(position, std::vector<glm::u8vec4>(), bounds);
    }
}

void makePKD(std::vector<glm::vec3>& position, std::vector<glm::u8vec4>& color, box3f const& bounds) {
    if (color.empty()) {
        recBuild(0, position.data(), position.size(), bounds);
    } else {
        recBuild(0, position.data(), position.size(), bounds, color.data());
    }
}

void makePKD(std::vector<glm::vec3>& position, std::vector<glm::u8vec4>& color, box3f const& bounds,
    unsigned int const begin, unsigned int const end) {
    if (color.empty()) {
        recBuild(0, position.data() + begin, end - begin, bounds);
    } else {
        recBuild(0, position.data() + begin, end - begin, bounds, color.data() + begin);
    }
}

std::vector<pkdlet> prePartition_inPlace(std::vector<glm::vec3>& particles, size_t maxSize, float radius, std::vector<glm::u8vec4>& color) {
    std::mutex resultMutex;
    std::vector<pkdlet> result;

    partitionRecursively(
        particles, 0ULL, particles.size(), [&](size_t begin, size_t end, bool force, box3f const& bounds) {
            /*bool makeLeaf() :*/
            const size_t size = end - begin;
            if (size > maxSize && !force)
                return false;

            pkdlet treelet;
            treelet.begin = begin;
            treelet.end = end;
            //treelet.bounds = extendBounds(particles, begin, end, radius);
            treelet.bounds = bounds;
            treelet.bounds.lower -= radius;
            treelet.bounds.upper += radius;

            std::lock_guard<std::mutex> lock(resultMutex);
            result.push_back(treelet);
            return true;
        }, color);

    return std::move(result);
}
} // namespace megamol::datatools
