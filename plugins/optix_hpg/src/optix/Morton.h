#pragma once

#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <glm/glm.hpp>

#include "box.h"
#include "particle.h"

namespace megamol::optix_hpg {
using MortonCode = glm::uvec3;

// https://stackoverflow.com/questions/49748864/morton-reverse-encoding-for-a-3d-grid

/* Morton encoding in binary (components 21-bit: 0..2097151)
                0zyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyx */
#define BITMASK_0000000001000001000001000001000001000001000001000001000001000001 UINT64_C(18300341342965825)
#define BITMASK_0000001000001000001000001000001000001000001000001000001000001000 UINT64_C(146402730743726600)
#define BITMASK_0001000000000000000000000000000000000000000000000000000000000000 UINT64_C(1152921504606846976)
/*              0000000ccc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc */
#define BITMASK_0000000000000011000000000011000000000011000000000011000000000011 UINT64_C(844631138906115)
#define BITMASK_0000000111000000000011000000000011000000000011000000000011000000 UINT64_C(126113986927919296)
/*              00000000000ccccc00000000cccc00000000cccc00000000cccc00000000cccc */
#define BITMASK_0000000000000000000000000000000000001111000000000000000000001111 UINT64_C(251658255)
#define BITMASK_0000000000000000000000001111000000000000000000001111000000000000 UINT64_C(1030792212480)
#define BITMASK_0000000000011111000000000000000000000000000000000000000000000000 UINT64_C(8725724278030336)
/*              000000000000000000000000000ccccccccccccc0000000000000000cccccccc */
#define BITMASK_0000000000000000000000000000000000000000000000000000000011111111 UINT64_C(255)
#define BITMASK_0000000000000000000000000001111111111111000000000000000000000000 UINT64_C(137422176256)
/*                                                         ccccccccccccccccccccc */
#define BITMASK_21BITS UINT64_C(2097151)


static inline void morton_decode(uint64_t m, uint32_t& xto, uint32_t& yto, uint32_t& zto) {
    constexpr uint64_t const mask0 = 0b0000000001000001000001000001000001000001000001000001000001000001,
                             mask1 = 0b0000001000001000001000001000001000001000001000001000001000001000,
                             mask2 = 0b0001000000000000000000000000000000000000000000000000000000000000,
                             mask3 = 0b0000000000000011000000000011000000000011000000000011000000000011,
                             mask4 = 0b0000000111000000000011000000000011000000000011000000000011000000,
                             mask5 = 0b0000000000000000000000000000000000001111000000000000000000001111,
                             mask6 = 0b0000000000000000000000001111000000000000000000001111000000000000,
                             mask7 = 0b0000000000011111000000000000000000000000000000000000000000000000,
                             mask8 = 0b0000000000000000000000000000000000000000000000000000000011111111,
                             mask9 = 0b0000000000000000000000000001111111111111000000000000000000000000;
    uint64_t x = m, y = m >> 1, z = m >> 2;

    /* 000c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c */
    x = (x & mask0) | ((x & mask1) >> 2) | ((x & mask2) >> 4);
    y = (y & mask0) | ((y & mask1) >> 2) | ((y & mask2) >> 4);
    z = (z & mask0) | ((z & mask1) >> 2) | ((z & mask2) >> 4);
    /* 0000000ccc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc */
    x = (x & mask3) | ((x & mask4) >> 4);
    y = (y & mask3) | ((y & mask4) >> 4);
    z = (z & mask3) | ((z & mask4) >> 4);
    /* 00000000000ccccc00000000cccc00000000cccc00000000cccc00000000cccc */
    x = (x & mask5) | ((x & mask6) >> 8) | ((x & mask7) >> 16);
    y = (y & mask5) | ((y & mask6) >> 8) | ((y & mask7) >> 16);
    z = (z & mask5) | ((z & mask6) >> 8) | ((z & mask7) >> 16);
    /* 000000000000000000000000000ccccccccccccc0000000000000000cccccccc */
    x = (x & mask8) | ((x & mask9) >> 16);
    y = (y & mask8) | ((y & mask9) >> 16);
    z = (z & mask8) | ((z & mask9) >> 16);
    /* 0000000000000000000000000000000000000000000ccccccccccccccccccccc */

    xto = x;

    yto = y;

    zto = z;
}


static inline uint64_t morton_encode(uint32_t xsrc, uint32_t ysrc, uint32_t zsrc) {
    constexpr uint64_t const mask0 = 0b0000000001000001000001000001000001000001000001000001000001000001,
                             mask1 = 0b0000001000001000001000001000001000001000001000001000001000001000,
                             mask2 = 0b0001000000000000000000000000000000000000000000000000000000000000,
                             mask3 = 0b0000000000000011000000000011000000000011000000000011000000000011,
                             mask4 = 0b0000000111000000000011000000000011000000000011000000000011000000,
                             mask5 = 0b0000000000000000000000000000000000001111000000000000000000001111,
                             mask6 = 0b0000000000000000000000001111000000000000000000001111000000000000,
                             mask7 = 0b0000000000011111000000000000000000000000000000000000000000000000,
                             mask8 = 0b0000000000000000000000000000000000000000000000000000000011111111,
                             mask9 = 0b0000000000000000000000000001111111111111000000000000000000000000;
    uint64_t x = xsrc, y = ysrc, z = zsrc;
    /* 0000000000000000000000000000000000000000000ccccccccccccccccccccc */
    x = (x & mask8) | ((x << 16) & mask9);
    y = (y & mask8) | ((y << 16) & mask9);
    z = (z & mask8) | ((z << 16) & mask9);
    /* 000000000000000000000000000ccccccccccccc0000000000000000cccccccc */
    x = (x & mask5) | ((x << 8) & mask6) | ((x << 16) & mask7);
    y = (y & mask5) | ((y << 8) & mask6) | ((y << 16) & mask7);
    z = (z & mask5) | ((z << 8) & mask6) | ((z << 16) & mask7);
    /* 00000000000ccccc00000000cccc00000000cccc00000000cccc00000000cccc */
    x = (x & mask3) | ((x << 4) & mask4);
    y = (y & mask3) | ((y << 4) & mask4);
    z = (z & mask3) | ((z << 4) & mask4);
    /* 0000000ccc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc */
    x = (x & mask0) | ((x << 2) & mask1) | ((x << 4) & mask2);
    y = (y & mask0) | ((y << 2) & mask1) | ((y << 4) & mask2);
    z = (z & mask0) | ((z << 2) & mask1) | ((z << 4) & mask2);
    /* 000c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c */
    return x | (y << 1) | (z << 2);
}

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
    std::vector<device::PKDParticle> const& data, device::box3f const& bounds) {
    std::vector<MortonCode> codes(data.size());

    constexpr uint64_t const factor = 1 << 21;
    //constexpr uint64_t const factor = 1 << 15;
    constexpr auto const dfactor = static_cast<double>(factor);

    auto const span = glm::dvec3(bounds.span());
    auto const lower = glm::dvec3(bounds.lower);

    for (size_t i = 0; i < data.size(); ++i) {
        auto const pos = (glm::dvec3(data[i].pos) - lower) / span;
        codes[i] = glm::uvec3(pos * dfactor);
    }

    std::vector<std::pair<uint64_t, uint64_t>> mc(codes.size());
    for (size_t i = 0; i < codes.size(); ++i) {
        mc[i] = std::make_pair(morton_encode(codes[i].x, codes[i].y, codes[i].z), i);
    }

    return mc;
}

std::vector<std::pair<uint64_t, uint64_t>> create_morton_codes(
    std::vector<device::PKDParticle> const& data, size_t begin, size_t end, device::box3f const& bounds) {
    auto const datasize = end - begin;
    std::vector<MortonCode> codes(datasize);

    constexpr uint64_t const factor = 1 << 21;
    //constexpr uint64_t const factor = 1 << 15;
    constexpr auto const dfactor = static_cast<double>(factor);

    auto const span = glm::dvec3(bounds.span());
    auto const lower = glm::dvec3(bounds.lower);

    for (size_t i = begin; i < end; ++i) {
        auto const pos = (glm::dvec3(data[i].pos) - lower) / span;
        codes[i - begin] = glm::uvec3(pos * dfactor);
    }

    std::vector<std::pair<uint64_t, uint64_t>> mc(datasize);
    for (size_t i = begin; i < end; ++i) {
        mc[i - begin] = std::make_pair(morton_encode(codes[i - begin].x, codes[i - begin].y, codes[i - begin].z), i);
    }

    return mc;
}

void sort_morton_codes(std::vector<std::pair<uint64_t, uint64_t>>& codes) {
    std::sort(codes.begin(), codes.end(), [](auto const& lhs, auto const& rhs) { return lhs.first < rhs.first; });
}

std::tuple<std::vector<std::pair<uint64_t, uint64_t>>, std::vector<uint64_t>, std::vector<device::PKDParticle>>
mask_morton_codes(
    std::vector<std::pair<uint64_t, uint64_t>> const& codes, std::vector<device::PKDParticle> const& data) {
    //constexpr uint32_t const mask = 0b11111111111111110000000000;
    //constexpr uint32_t const mask = 0b1111111111;
    //constexpr uint64_t const mask = 0b111111111111111111000000000000000000000000000000000000000000000;
    constexpr uint64_t const mask = 0b111111111111111000000000000000000000000000000000000000000000000;

    //constexpr uint64_t const mask = 0b111111111111111000000000000000000000000000000;


    std::unordered_map<uint64_t, uint64_t> grid;

    for (size_t i = 0; i < codes.size(); ++i) {
        ++grid[(codes[i].first & mask) >> 48];
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
        auto const idx = --grid[(codes[i].first & mask) >> 48];
        sorted_codes[idx] = codes[i].first;
        sorted_data[idx] = data[codes[i].second];
    }

    return std::make_tuple(cells, sorted_codes, sorted_data);
}

void convert_morton_treelet(
    device::PKDlet const& treelet, std::vector<device::PKDParticle> const& data, device::box3f const& global_bounds) {
    auto const codes = create_morton_codes(data, treelet.begin, treelet.end, global_bounds);

    constexpr uint64_t const factor = 1 << 21;
    constexpr uint64_t const mask = 0b111111111111111000000000000000000000000000000000000000000000000;
    constexpr uint64_t const mask2 = 0b000000000000000111111111111111111111111111111000000000000000000;

    std::unordered_map<uint64_t, uint64_t> grid;

    for (size_t i = 0; i < codes.size(); ++i) {
        ++grid[(codes[i].first & mask) >> 48];
    }

    /*uint32_t x, y, z;
    morton_decode(grid.begin()->first, x, y, z);

    glm::vec3 basePos(
        x / static_cast<double>(factor), y / static_cast<double>(factor), z / static_cast<double>(factor));
    basePos *= global_bounds.span();
    basePos += global_bounds.lower;*/
    std::vector<glm::vec3> recon_data(codes.size());
    for (size_t i = 0; i < codes.size(); ++i) {
        auto const code = (codes[i].first >> 18) << 18;
        uint32_t x, y, z;
        morton_decode(code, x, y, z);
        glm::vec3 basePos(
            x / static_cast<double>(factor), y / static_cast<double>(factor), z / static_cast<double>(factor));
        basePos *= global_bounds.span();
        basePos += global_bounds.lower;
        recon_data[i] = basePos;
    }

    std::vector<glm::dvec3> diffs(codes.size());
    for (size_t i = treelet.begin; i < treelet.end; ++i) {
        diffs[i - treelet.begin] = glm::abs(glm::dvec3(data[i].pos) - glm::dvec3(recon_data[i - treelet.begin]));
    }
}
} // namespace megamol::optix_hpg
