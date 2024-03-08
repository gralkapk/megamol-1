#pragma once

#include <glm/glm.hpp>

namespace megamol::optix_hpg {
template<typename FType, typename BaseType, BaseType D>
class FixedPoint {
    static_assert(std::is_floating_point_v<FType>, "FType must be floating point type");
    static_assert(std::is_integral_v<BaseType>, "BaseType must be integral type");

public:
    constexpr static BaseType factor = 1 << D;
    constexpr static BaseType half_factor = 1 << (D - 1);

    FixedPoint() = default;

    explicit FixedPoint(FType val) {
        *this = val;
    }

    FixedPoint(FixedPoint const& rhs) = default;
    FixedPoint& operator=(FixedPoint const& rhs) = default;
    FixedPoint(FixedPoint&& rhs) = default;
    FixedPoint& operator=(FixedPoint&& rhs) = default;

    /*FixedPoint(FixedPoint const& rhs) {
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
        return *this;
    }*/

    FixedPoint& operator=(FType val) {
        value_ = static_cast<BaseType>(val * factor);
        return *this;
    }

    FixedPoint operator+(FixedPoint const& rhs) const {
        FixedPoint ret;
        ret.value_ = value_ + rhs.value_;
        return ret;
    }

    FixedPoint operator-(FixedPoint const& rhs) const {
        FixedPoint ret;
        ret.value_ = value_ - rhs.value_;
        return ret;
    }

    FixedPoint operator*(FixedPoint const& rhs) const {
        FixedPoint ret;
        ret.value_ = ((value_ * rhs.value_) + half_factor) >> D;
        return ret;
    }

    FixedPoint operator/(FixedPoint const& rhs) const {
        FixedPoint ret;
        ret.value_ = (value_ << D) / rhs.value_;
        return ret;
    }

    operator FType() const {
        return static_cast<FType>(value_) / static_cast<FType>(factor);
    }

    operator BaseType() const {
        return value_;
    }

private:
    BaseType value_;
};

using decvec3 = glm::vec<3, FixedPoint<float, unsigned, 8>>;
} // namespace megamol::optix_hpg
