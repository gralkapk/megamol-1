#pragma once

#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <glm/glm.hpp>

#include "box.h"
#include "morton_util.h"
#include "particle.h"

namespace megamol::optix_hpg {
using MortonCode = glm::uvec3;


// wikipedia

bool less_msb(unsigned int lhs, unsigned int rhs) {
    return lhs < rhs && lhs < (lhs ^ rhs);
}

bool cmp_mc(MortonCode const& lhs, MortonCode const& rhs) {
    int msd = 0;
    for (int i = 1; i < 3; ++i) {
        if (less_msb(lhs[msd] ^ rhs[msd], lhs[i] ^ rhs[i])) {
            msd = i;
        }
    }
    return lhs[msd] < rhs[msd];
}

std::vector<std::pair<uint64_t, uint64_t>> create_morton_codes(
    std::vector<device::PKDParticle> const& data, device::box3f const& bounds, device::MortonConfig const& config) {
    std::vector<MortonCode> codes(data.size());

    //constexpr uint64_t const factor = 1 << 21;
    //constexpr uint64_t const factor = 1 << 15;
    auto const dfactor = static_cast<double>(config.factor);

    auto const span = glm::dvec3(bounds.span());
    auto const lower = glm::dvec3(bounds.lower);

    for (size_t i = 0; i < data.size(); ++i) {
        auto const pos = (glm::dvec3(data[i].pos) - lower) / span;
        codes[i] = glm::uvec3(pos * dfactor);
    }

    std::vector<std::pair<uint64_t, uint64_t>> mc(codes.size());
    for (size_t i = 0; i < codes.size(); ++i) {
        mc[i] = std::make_pair(device::morton_encode(codes[i].x, codes[i].y, codes[i].z), i);
    }

    return mc;
}

std::vector<std::pair<uint64_t, uint64_t>> create_morton_codes(std::vector<device::PKDParticle> const& data,
    size_t begin, size_t end, device::box3f const& bounds, device::MortonConfig const& config) {
    auto const datasize = end - begin;
    std::vector<MortonCode> codes(datasize);

    //constexpr uint64_t const factor = 1 << 21;
    //constexpr uint64_t const factor = 1 << 15;
    auto const dfactor = static_cast<double>(config.factor);

    auto const span = glm::dvec3(bounds.span());
    auto const lower = glm::dvec3(bounds.lower);

    for (size_t i = begin; i < end; ++i) {
        auto const pos = (glm::dvec3(data[i].pos) - lower) / span;
        codes[i - begin] = glm::uvec3(pos * dfactor);
    }

    std::vector<std::pair<uint64_t, uint64_t>> mc(datasize);
    for (size_t i = begin; i < end; ++i) {
        mc[i - begin] =
            std::make_pair(device::morton_encode(codes[i - begin].x, codes[i - begin].y, codes[i - begin].z), i);
    }

    return mc;
}

void sort_morton_codes(std::vector<std::pair<uint64_t, uint64_t>>& codes) {
    std::sort(codes.begin(), codes.end(), [](auto const& lhs, auto const& rhs) { return lhs.first < rhs.first; });
}

std::tuple<std::vector<std::pair<uint64_t, uint64_t>>, std::vector<uint64_t>, std::vector<device::PKDParticle>>
mask_morton_codes(std::vector<std::pair<uint64_t, uint64_t>> const& codes, std::vector<device::PKDParticle> const& data,
    device::MortonConfig const& config) {
    //constexpr uint32_t const mask = 0b11111111111111110000000000;
    //constexpr uint32_t const mask = 0b1111111111;
    //constexpr uint64_t const mask = 0b111111111111111111000000000000000000000000000000000000000000000;
    //constexpr uint64_t const mask = 0b111111111111111000000000000000000000000000000000000000000000000;

    //constexpr uint64_t const mask = 0b111111111111111000000000000000000000000000000;


    std::unordered_map<uint64_t, uint64_t> grid;

    for (size_t i = 0; i < codes.size(); ++i) {
        ++grid[(codes[i].first & config.mask) >> config.offset];
        //grid[i] = (codes[i] & mask) >> 30;
    }

    //std::sort(grid.begin(), grid.end(), [](auto const& lhs, auto const& rhs) { return lhs.x < rhs.x; });

    //grid.erase(std::unique(grid.begin(), grid.end()), grid.end());

    std::vector<std::pair<uint64_t, uint64_t>> grid_helper(grid.begin(), grid.end());
    std::sort(
        grid_helper.begin(), grid_helper.end(), [](auto const& lhs, auto const& rhs) { return lhs.first < rhs.first; });

    std::vector<std::pair<uint64_t, uint64_t>> cells(grid_helper.size());
    cells[0].second = grid_helper[0].second;
    for (size_t i = 1; i < grid_helper.size(); ++i) {
        grid_helper[i].second += grid_helper[i - 1].second;
        cells[i].first = grid_helper[i - 1].second;
        cells[i].second = grid_helper[i].second;
    }

    grid.clear();
    grid.reserve(grid_helper.size());
    grid.insert(grid_helper.begin(), grid_helper.end());

    std::vector<uint64_t> sorted_codes(codes.size());
    std::vector<device::PKDParticle> sorted_data(data.size());

    for (size_t i = 0; i < codes.size(); ++i) {
        auto const idx = --grid[(codes[i].first & config.mask) >> config.offset];
        sorted_codes[idx] = codes[i].first;
        sorted_data[idx] = data[codes[i].second];
    }

    return std::make_tuple(cells, sorted_codes, sorted_data);
}

void convert_morton_treelet(device::PKDlet const& treelet, std::vector<device::PKDParticle> const& data,
    device::C2PKDlet& ctreelet, std::vector<device::C2PKDParticle>& cparticles, device::box3f const& global_bounds,
    device::MortonConfig const& config) {
    auto const codes = create_morton_codes(data, treelet.begin, treelet.end, global_bounds, config);

    //constexpr uint64_t const factor = 1 << 21;
    //constexpr uint64_t const mask = 0b111111111111111000000000000000000000000000000000000000000000000;
    //constexpr uint64_t const mask2 = 0b000000000000000111111111111111111111111111111000000000000000000;

    std::unordered_map<uint64_t, uint64_t> grid;

    for (size_t i = 0; i < codes.size(); ++i) {
        ++grid[(codes[i].first & config.mask) >> config.offset];
    }

    auto const global_prefix = (codes[0].first & config.mask) >> config.offset;
    ctreelet.prefix = global_prefix;
    ctreelet.begin = treelet.begin;
    ctreelet.end = treelet.end;
    ctreelet.bounds = treelet.bounds;

    auto const span = global_bounds.span();
    auto const lower = global_bounds.lower;

    std::vector<glm::vec3> recon_data(codes.size());
    for (size_t i = 0; i < codes.size(); ++i) {
        auto const code = codes[i].first >> config.code_offset;
        auto const prefix = (codes[i].first & config.mask) >> config.offset;
        if (global_prefix != prefix) {
            throw std::runtime_error("unexpected prefix");
        }

        cparticles[i + treelet.begin].code = code;

        /*auto const combined_code = (static_cast<uint64_t>(cparticles[i + treelet.begin].code) << config.code_offset) +
                                   (static_cast<uint64_t>(ctreelet.prefix) << config.offset);

        uint32_t x, y, z;
        device::morton_decode(combined_code, x, y, z);
        glm::vec3 basePos(x / static_cast<double>(config.factor), y / static_cast<double>(config.factor),
            z / static_cast<double>(config.factor));
        basePos *= global_bounds.span();
        basePos += global_bounds.lower;*/

        auto const basePos = cparticles[i + treelet.begin].from(
            ctreelet.prefix, span, lower, config.code_offset, config.offset, config.factor);

        recon_data[i] = basePos;
    }

    std::vector<glm::dvec3> diffs(codes.size());
    for (size_t i = treelet.begin; i < treelet.end; ++i) {
        diffs[i - treelet.begin] = glm::abs(glm::dvec3(data[i].pos) - glm::dvec3(recon_data[i - treelet.begin]));
    }
}

void adapt_morton_bbox(std::vector<device::C2PKDParticle> const& cparticles, device::C2PKDlet& treelet,
    device::box3f const& global_bounds, float const radius, device::MortonConfig const& config) {
    device::box3f bounds;

    auto const span = global_bounds.span();
    auto const lower = global_bounds.lower;

    for (unsigned int i = treelet.begin; i < treelet.end; ++i) {
        bounds.extend(
            cparticles[i].from(treelet.prefix, span, lower, config.code_offset, config.offset, config.factor));
    }

    bounds.lower -= radius;
    bounds.upper += radius;

    treelet.bounds = bounds;
}
} // namespace megamol::optix_hpg
