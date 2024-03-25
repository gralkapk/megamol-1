#pragma once

#include <cstddef>
#include <queue>
#include <tuple>
#include <vector>
#include <unordered_set>

#include "particle.h"

#include <cuda.h>

namespace megamol::optix_hpg {
inline bool tc_comp(std::tuple<int, int> const& lhs, std::tuple<int, int> const& rhs) {
    return std::get<0>(lhs) > std::get<0>(rhs);
}

class TreeletCache {
public:
    TreeletCache(std::size_t cache_size, std::vector<device::PKDlet> const& treelets);

    ~TreeletCache();

    CUdeviceptr Alloc(int treeletID);

    void Free(int treeletID);

private:
    std::priority_queue<std::tuple<int /*frame id*/, int /*treelet id*/>, std::vector<std::tuple<int, int>>,
        decltype(tc_comp)>
        lru;

    CUdeviceptr vptr_ = 0;

    std::size_t vptr_res_size_ = 0;

    std::vector<std::size_t> treelet_ptr_offs_;

    std::vector<std::size_t> treelet_ptr_szs_;

    std::size_t gran_ = 0;

    std::vector<CUdeviceptr> treelet_ptrs_;

    std::unordered_set<int> committed_treelets_;
};
} // namespace megamol::optix_hpg
