#pragma once

#include <glm/glm.hpp>

namespace megamol::optix_hpg {
template<typename FType, typename BaseType, BaseType D>
class FixedPoint {
    static_assert(std::is_floating_point_v<FType>, "FType must be floating point type");
    static_assert(std::is_integral_v<BaseType>, "BaseType must be integral type");

public:
    constexpr static BaseType factor = 1 << D;

    explicit FixedPoint(FType val) {
        *this = val;
    }

    FixedPoint(FixedPoint const& rhs) {
        *this = rhs;
    }

    FixedPoint& operator=(FixedPoint const& rhs) {
        value_ = rhs.value_;
    }

    FixedPoint(FixedPoint&& rhs) {
        *this = std::move(rhs);
    }

    FixedPoint& operator=(FixedPoint&& rhs) {
        value_ = std::exchange(rhs.value_, static_cast<BaseType>(0));
    }

    FixedPoint& operator=(FType val) {
        value_ = static_cast<BaseType>(val * factor);
        return *this;
    }

    FixedPoint operator+(FixedPoint const& rhs) const {
        return value_ + rhs.value_;
    }

    FixedPoint operator-(FixedPoint const& rhs) const {
        return value_ - rhs.value_;
    }

    FixedPoint operator*(FixedPoint const& rhs) const {
        return (value_ * rhs.value_) >> D;
    }

    FixedPoint operator/(FixedPoint const& rhs) const {
        return (value_ / rhs.value_) << D;
    }

    operator FType() const {
        return static_cast<FType>(value_) / static_cast<FType>(factor);
    }

private:
    BaseType value_;
};

using decvec3 = glm::vec<3, FixedPoint<float, unsigned, 8>>;
} // namespace megamol::optix_hpg
