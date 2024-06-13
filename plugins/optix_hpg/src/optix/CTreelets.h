#pragma once

#include <unordered_map>
#include <vector>

#include "FTreelets.h"
#include "PKDUtils.h"

#include "particle.h"

#include "tbb/parallel_for.h"

#include "omp.h"

namespace megamol::optix_hpg {
inline void norm_at_bounds(std::vector<device::PKDParticle>& data, device::box3f const& bounds) {
    auto const span = glm::dvec3(bounds.span());
    auto const lower = glm::dvec3(bounds.lower);

    tbb::parallel_for(
        (size_t) 0, data.size(), [&](size_t p_id) { data[p_id].pos = (glm::dvec3(data[p_id].pos) - lower) / span; });
}

inline void ctreelets_partition(
    std::vector<device::PKDParticle>& data, device::box3f const& bounds, float radius, size_t maxSize) {
    auto fparticles = convert_to_fparticles(data, bounds);

    auto const max_threads = omp_get_max_threads();
    std::vector<device::box3u32> thread_box(max_threads);
#pragma omp parallel for
    for (int64_t i = 0; i < fparticles.size(); ++i) {
        thread_box[omp_get_thread_num()].extend(fparticles[i].pos);
    }

    for (int i = 1; i < thread_box.size(); ++i) {
        thread_box[0].extend(thread_box[i]);
    }

    auto const& fbounds = thread_box[0];

    std::mutex resultMutex;
    std::vector<device::FPKDLet> treelets;
    std::vector<glm::uvec3> prefix_counter;

    constexpr uint64_t const factor = 1 << 23;
    auto const span = glm::dvec3(bounds.span());
    glm::uvec3 scaled_radius =
        glm::uvec3((radius / span.x) * factor, (radius / span.y) * factor, (radius / span.z) * factor);

    partitionRecursively<device::FPKDParticle, device::box3u32>(
        fparticles, 0ULL, fparticles.size(), [&](size_t begin, size_t end, bool force, device::box3u32 const& bounds) {
            /*bool makeLeaf() :*/
            const size_t size = end - begin;
            if (size > maxSize && !force)
                return false;

            std::unordered_map<uint32_t, uint32_t> x_counter;
            std::unordered_map<uint32_t, uint32_t> y_counter;
            std::unordered_map<uint32_t, uint32_t> z_counter;
            for (size_t i = begin; i < end; ++i) {
                ++x_counter[(fparticles[i].pos.x - bounds.lower.x) >> 12];
                ++y_counter[(fparticles[i].pos.y - bounds.lower.y) >> 12];
                ++z_counter[(fparticles[i].pos.z - bounds.lower.z) >> 12];
            }
            /*if (x_counter.size() > 1 || y_counter.size() > 1 || z_counter.size() > 1)
                return false;*/

            device::FPKDLet treelet;
            treelet.begin = begin;
            treelet.end = end;
            treelet.bounds = bounds;
            /*treelet.bounds.lower -= scaled_radius;
            treelet.bounds.upper += scaled_radius;*/

            std::lock_guard<std::mutex> lock(resultMutex);
            treelets.push_back(treelet);
            prefix_counter.push_back(glm::uvec3(x_counter.size(), y_counter.size(), z_counter.size()));
            return true;
        });

    auto const pc = std::count_if(
        prefix_counter.begin(), prefix_counter.end(), [](auto const& el) { return el.x > 1 || el.y > 1 || el.z > 1; });
}
} // namespace megamol::optix_hpg
