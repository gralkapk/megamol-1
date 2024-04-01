#include "TreeletCache.h"

#include "Utils.h"

namespace megamol::optix_hpg {
TreeletCache::TreeletCache(std::size_t cache_size, std::vector<device::PKDlet> const& treelets) {
    // get the required granularity for the chunks
    CUmemAllocationProp prop = {};
    prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
    prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    prop.location.id = 0;
    CUDA_CHECK_ERROR(cuMemGetAllocationGranularity(&gran_, &prop, CU_MEM_ALLOC_GRANULARITY_RECOMMENDED));

    for (int tID = 0; tID < treelets.size(); ++tID) {
        auto const& treelet = treelets[tID];
        auto const byte_begin = treelet.begin * sizeof(device::QPKDParticle);
        auto const byte_end = treelet.end * sizeof(device::QPKDParticle);
        auto const begin_chunk_idx = byte_begin / gran_;
        auto const end_chunk_idx = byte_end / gran_;
        for (auto ci = begin_chunk_idx; ci <= end_chunk_idx; ++ci) {
            treelets_to_chunks_[tID].push_back(ci);
            chunks_to_treelets_[ci].push_back(tID);
        }
        treelet_ptr_offs_[tID] = byte_begin;
    }

    vptr_res_size_ = chunks_to_treelets_.size() * gran_;
    CUDA_CHECK_ERROR(cuMemAddressReserve(&vptr_, vptr_res_size_, 0, 0, 0));

    #if 0
    treelet_ptr_offs_.reserve(treelets.size());
    treelet_ptr_szs_.reserve(treelets.size());
    treelet_ptrs_.resize(treelets.size());
    CUmemAllocationProp prop = {};
    prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
    prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    prop.location.id = 0;
    CUDA_CHECK_ERROR(cuMemGetAllocationGranularity(&gran_, &prop, CU_MEM_ALLOC_GRANULARITY_RECOMMENDED));
    std::size_t padded_total_treelet_size = 0;
    for (auto const& t : treelets) {
        auto const factor = (t.end - t.begin) * sizeof(device::QPKDParticle) / static_cast<double>(gran_);
        treelet_ptr_offs_.push_back(padded_total_treelet_size);
        padded_total_treelet_size += std::ceil(factor) * gran_;
        treelet_ptr_szs_.push_back(std::ceil(factor) * gran_);
    }
    auto const factor = padded_total_treelet_size / static_cast<double>(gran_);
    vptr_res_size_ = std::ceil(factor) * gran_;
    CUDA_CHECK_ERROR(cuMemAddressReserve(&vptr_, vptr_res_size_, 0, 0, 0));
    #endif
}

TreeletCache::~TreeletCache() {
    CUDA_CHECK_ERROR(cuMemAddressFree(vptr_, vptr_res_size_));
}

CUdeviceptr TreeletCache::Alloc(int treeletID) {
    if (committed_treelets_.find(treeletID) != committed_treelets_.end()) {
        return 0;
    }

    auto const& req_chunks = treelets_to_chunks_.at(treeletID);
    for (auto const& rc : req_chunks) {
        if (committed_chunks_.find(rc) != committed_chunks_.end())
            continue;
        CUdeviceptr handle = 0;
        CUmemAllocationProp prop = {};
        prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
        prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
        prop.location.id = 0;
        CUDA_CHECK_ERROR(cuMemCreate(&handle, gran_, &prop, 0));
        CUDA_CHECK_ERROR(cuMemMap(vptr_ + rc * gran_, gran_, 0, handle, 0));
        CUmemAccessDesc accessDesc = {};
        accessDesc.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
        accessDesc.location.id = 0;
        accessDesc.flags = CU_MEM_ACCESS_FLAGS_PROT_READWRITE;
        CUDA_CHECK_ERROR(cuMemSetAccess(vptr_ + rc * gran_, gran_, &accessDesc, 1));
        chunk_ptrs_[rc] = handle;
        committed_chunks_.insert(rc);
    }
    committed_treelets_.insert(treeletID);
    return vptr_ + treelet_ptr_offs_[treeletID];

#if 0
    CUdeviceptr handle = 0;
    CUmemAllocationProp prop = {};
    prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
    prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    prop.location.id = 0;
    CUDA_CHECK_ERROR(cuMemCreate(&handle, treelet_ptr_szs_[treeletID], &prop, 0));
    CUDA_CHECK_ERROR(cuMemMap(vptr_ + treelet_ptr_offs_[treeletID], treelet_ptr_szs_[treeletID], 0, handle, 0));
    CUmemAccessDesc accessDesc = {};
    accessDesc.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    accessDesc.location.id = 0;
    accessDesc.flags = CU_MEM_ACCESS_FLAGS_PROT_READWRITE;
    CUDA_CHECK_ERROR(cuMemSetAccess(vptr_ + treelet_ptr_offs_[treeletID], treelet_ptr_szs_[treeletID], &accessDesc, 1));
    treelet_ptrs_[treeletID] = handle;
    committed_treelets_.insert(treeletID);
    return vptr_ + treelet_ptr_offs_[treeletID];
#endif
}

void TreeletCache::Free(int treeletID) {
    CUDA_CHECK_ERROR(cuMemUnmap(vptr_ + treelet_ptr_offs_[treeletID], treelet_ptr_szs_[treeletID]));
    CUDA_CHECK_ERROR(cuMemFree(treelet_ptrs_[treeletID]));
    committed_treelets_.erase(treeletID);
}
} // namespace megamol::optix_hpg
