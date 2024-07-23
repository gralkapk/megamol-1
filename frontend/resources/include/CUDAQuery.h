#pragma once

#include <AnyQuery.h>

namespace megamol::frontend_resources::performance {
///<summary>
/// Wrapper for CUDA timer query.
/// </summary>
class CUDAQuery : public AnyQuery {
public:
    CUDAQuery();

    ~CUDAQuery() override;

    /// <summary>
    /// Set timestamp query.
    /// </summary>
    void Counter(void* userData = nullptr) override;

    std::shared_ptr<AnyQuery> MakeAnother() override {
        return std::make_shared<CUDAQuery>();
    }

    /// <summary>
    /// Try to retrieve the timestamp.
    /// Does not wait if value is not ready.
    /// After successful retrieval will return acquired timestamp and not try again.
    /// </summary>
    /// <returns>Queried timestamp or zero if value is not ready</returns>
    time_point GetNW() override;
};
} // namespace megamol::frontend_resources::performance
