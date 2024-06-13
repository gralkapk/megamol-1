#pragma once

#include <iostream>
#include <stack>
#include <tuple>
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

    /*constexpr uint64_t const factor = 1 << 23;
    auto const span = glm::dvec3(bounds.span());
    glm::uvec3 scaled_radius =
        glm::uvec3((radius / span.x) * factor, (radius / span.y) * factor, (radius / span.z) * factor);*/

    std::stack<std::tuple<size_t, size_t, int, int, int>> remaining_jobs;
    int offset = 10;

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
                ++x_counter[(fparticles[i].pos.x - bounds.lower.x) >> offset];
                ++y_counter[(fparticles[i].pos.y - bounds.lower.y) >> offset];
                ++z_counter[(fparticles[i].pos.z - bounds.lower.z) >> offset];
            }
            /*if (x_counter.size() > 1 || y_counter.size() > 1 || z_counter.size() > 1)
                return false;*/
            if (x_counter.size() <= 1 && y_counter.size() <= 1 && z_counter.size() <= 1) {
                device::FPKDLet treelet;
                treelet.begin = begin;
                treelet.end = end;
                treelet.bounds = bounds;
                /*treelet.bounds.lower -= scaled_radius;
                treelet.bounds.upper += scaled_radius;*/
                treelet.offset[0] = offset;
                treelet.offset[1] = offset;
                treelet.offset[2] = offset;

                std::lock_guard<std::mutex> lock(resultMutex);
                treelets.push_back(treelet);
                prefix_counter.push_back(glm::uvec3(x_counter.size(), y_counter.size(), z_counter.size()));
            } else {
                auto offset_x = offset;
                auto offset_y = offset;
                auto offset_z = offset;
                if (x_counter.size() > 1) {
                    ++offset_x;
                }
                if (y_counter.size() > 1) {
                    ++offset_y;
                }
                if (z_counter.size() > 1) {
                    ++offset_z;
                }

                std::lock_guard<std::mutex> lock(resultMutex);
                remaining_jobs.push(std::make_tuple(begin, end, offset_x, offset_y, offset_z));
            }
            return true;
        });

    while (!remaining_jobs.empty()) {
        auto rj = remaining_jobs.top();
        remaining_jobs.pop();
        partitionRecursively<device::FPKDParticle, device::box3u32>(fparticles, std::get<0>(rj), std::get<1>(rj),
            [&](size_t begin, size_t end, bool force, device::box3u32 const& bounds) {
                /*bool makeLeaf() :*/
                const size_t size = end - begin;
                if (size > maxSize && !force)
                    return false;

                auto const rj_offset_x = std::get<2>(rj);
                auto const rj_offset_y = std::get<3>(rj);
                auto const rj_offset_z = std::get<4>(rj);
                std::unordered_map<uint32_t, uint32_t> x_counter;
                std::unordered_map<uint32_t, uint32_t> y_counter;
                std::unordered_map<uint32_t, uint32_t> z_counter;
                for (size_t i = begin; i < end; ++i) {
                    ++x_counter[(fparticles[i].pos.x - bounds.lower.x) >> rj_offset_x];
                    ++y_counter[(fparticles[i].pos.y - bounds.lower.y) >> rj_offset_y];
                    ++z_counter[(fparticles[i].pos.z - bounds.lower.z) >> rj_offset_z];
                }
                /*if (x_counter.size() > 1 || y_counter.size() > 1 || z_counter.size() > 1)
                    return false;*/
                if (x_counter.size() <= 1 && y_counter.size() <= 1 && z_counter.size() <= 1) {
                    device::FPKDLet treelet;
                    treelet.begin = begin;
                    treelet.end = end;
                    treelet.bounds = bounds;
                    /*treelet.bounds.lower -= scaled_radius;
                    treelet.bounds.upper += scaled_radius;*/
                    treelet.offset[0] = rj_offset_x;
                    treelet.offset[1] = rj_offset_y;
                    treelet.offset[2] = rj_offset_z;

                    std::lock_guard<std::mutex> lock(resultMutex);
                    treelets.push_back(treelet);
                    prefix_counter.push_back(glm::uvec3(x_counter.size(), y_counter.size(), z_counter.size()));
                } else {
                    auto offset_x = rj_offset_x;
                    auto offset_y = rj_offset_y;
                    auto offset_z = rj_offset_z;
                    if (x_counter.size() > 1) {
                        ++offset_x;
                    }
                    if (y_counter.size() > 1) {
                        ++offset_y;
                    }
                    if (z_counter.size() > 1) {
                        ++offset_z;
                    }
                    std::lock_guard<std::mutex> lock(resultMutex);
                    remaining_jobs.push(std::make_tuple(begin, end, offset_x, offset_y, offset_z));
                }
                return true;
            }, true);
    }

    auto const pc = std::count_if(
        prefix_counter.begin(), prefix_counter.end(), [](auto const& el) { return el.x > 1 || el.y > 1 || el.z > 1; });
    auto const toc_x =
        std::count_if(treelets.begin(), treelets.end(), [&](auto const& el) { return el.offset[0] > offset; });
    auto const toc_y =
        std::count_if(treelets.begin(), treelets.end(), [&](auto const& el) { return el.offset[1] > offset; });
    auto const toc_z =
        std::count_if(treelets.begin(), treelets.end(), [&](auto const& el) { return el.offset[2] > offset; });
    auto const max_to_x = std::max_element(treelets.begin(), treelets.end(),
        [](auto const& lhs, auto const& rhs) { return lhs.offset[0] < rhs.offset[0]; });
    auto const max_to_y = std::max_element(treelets.begin(), treelets.end(),
        [](auto const& lhs, auto const& rhs) { return lhs.offset[1] < rhs.offset[1]; });
    auto const max_to_z = std::max_element(treelets.begin(), treelets.end(),
        [](auto const& lhs, auto const& rhs) { return lhs.offset[2] < rhs.offset[2]; });

    std::cout << "pc " << pc << " toc_x " << toc_x << " toc_y " << toc_y << " toc_z " << toc_z << " max_to_x "
              << max_to_x->offset[0] << " max_to_y " << max_to_y->offset[1] << " max_to_z " << max_to_z->offset[2] << std::endl;
}
} // namespace megamol::optix_hpg
