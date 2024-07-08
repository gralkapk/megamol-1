#pragma once

#include <glm/glm.hpp>

#include "optix/utils_host.h"

namespace megamol {
namespace optix_hpg {

enum class QTreeletType { E5M15, E4M16, E5M15D, E4M16D };

namespace device {
struct QPKDlet {
    box3f bounds;

    size_t begin, end;

    glm::vec3 basePos;
};

struct QTParticle {
    unsigned char dim;
    char exp_x;
    char exp_y;
    char exp_z;
    unsigned char sign_x : 1;
    unsigned char sign_y : 1;
    unsigned char sign_z : 1;
    unsigned short m_x;
    unsigned short m_y;
    unsigned short m_z;
};

struct QTParticle_e5m15 {
    static constexpr int offset = 8;
    static constexpr bool dep = false;
    static constexpr char const* type_name = "e5m15";
    static constexpr int exp = 5;
    static constexpr bool e_sign = false;
    unsigned short exp_x : exp;
    unsigned short exp_y : exp;
    unsigned short exp_z : exp;
    unsigned short m_x : 15;
    unsigned short dim_x : 1;
    unsigned short m_y : 15;
    unsigned short dim_y : 1;
    unsigned short m_z : 15;
    unsigned short dim_z : 1;
    QTParticle getParticle() const {
        QTParticle qpart;
        qpart.exp_x = this->exp_x;
        qpart.exp_y = this->exp_y;
        qpart.exp_z = this->exp_z;
        qpart.m_x = this->m_x;
        qpart.m_y = this->m_y;
        qpart.m_z = this->m_z;
        if (this->dim_x) {
            qpart.dim = 0;
        }
        if (this->dim_y) {
            qpart.dim = 1;
        }
        if (this->dim_z) {
            qpart.dim = 2;
        }
        qpart.sign_x = 0;
        qpart.sign_y = 0;
        qpart.sign_z = 0;
        return qpart;
    }
};

struct QTParticle_e4m16 {
    static constexpr int offset = 7;
    static constexpr bool dep = false;
    static constexpr char const* type_name = "e4m16";
    static constexpr int exp = 4;
    static constexpr bool e_sign = false;
    unsigned short exp_x : exp;
    unsigned short dim_x : 1;
    unsigned short exp_y : exp;
    unsigned short dim_y : 1;
    unsigned short exp_z : exp;
    unsigned short dim_z : 1;
    unsigned short m_x : 16;
    unsigned short m_y : 16;
    unsigned short m_z : 16;
    QTParticle getParticle() const {
        QTParticle qpart;
        qpart.exp_x = this->exp_x;
        qpart.exp_y = this->exp_y;
        qpart.exp_z = this->exp_z;
        qpart.m_x = this->m_x;
        qpart.m_y = this->m_y;
        qpart.m_z = this->m_z;
        if (this->dim_x) {
            qpart.dim = 0;
        }
        if (this->dim_y) {
            qpart.dim = 1;
        }
        if (this->dim_z) {
            qpart.dim = 2;
        }
        qpart.sign_x = 0;
        qpart.sign_y = 0;
        qpart.sign_z = 0;
        return qpart;
    }
};

struct QTParticle_e5m15d {
    static constexpr int offset = 8;
    static constexpr bool dep = true;
    static constexpr char const* type_name = "e5m15d";
    static constexpr int exp = 5;
    static constexpr bool e_sign = false;
    unsigned short exp_x : exp;
    unsigned short exp_y : exp;
    unsigned short exp_z : exp;
    unsigned short dim_a : 1;
    unsigned short m_x : 15;
    unsigned short dim_b : 1;
    unsigned short m_y : 15;
    unsigned short sign_a : 1;
    unsigned short m_z : 15;
    unsigned short sign_b : 1;
    QTParticle getParticle(bool left_child, int sep_dim) const {
        QTParticle qpart;
        qpart.exp_x = this->exp_x;
        qpart.exp_y = this->exp_y;
        qpart.exp_z = this->exp_z;
        qpart.m_x = this->m_x;
        qpart.m_y = this->m_y;
        qpart.m_z = this->m_z;
        qpart.dim = (this->dim_b << 1) + this->dim_a;
        qpart.sign_x = sep_dim == 0 ? left_child : this->sign_a;
        qpart.sign_y = sep_dim == 1 ? left_child : (sep_dim == 0 ? this->sign_a : this->sign_b);
        qpart.sign_z = sep_dim == 2 ? left_child : this->sign_b;
        return qpart;
    }
};

struct QTParticle_e4m16d {
    static constexpr int offset = 7;
    static constexpr bool dep = true;
    static constexpr char const* type_name = "e4m16d";
    static constexpr int exp = 4;
    static constexpr bool e_sign = false;
    unsigned short exp_x : exp;
    unsigned short exp_y : exp;
    unsigned short exp_z : exp;
    unsigned short dim_a : 1;
    unsigned short dim_b : 1;
    unsigned short sign_a : 1;
    unsigned short sign_b : 1;
    unsigned short m_x : 16;
    unsigned short m_y : 16;
    unsigned short m_z : 16;
    QTParticle getParticle(bool left_child, int sep_dim) const {
        QTParticle qpart;
        qpart.exp_x = this->exp_x;
        qpart.exp_y = this->exp_y;
        qpart.exp_z = this->exp_z;
        qpart.m_x = this->m_x;
        qpart.m_y = this->m_y;
        qpart.m_z = this->m_z;
        qpart.dim = (this->dim_b << 1) + this->dim_a;
        qpart.sign_x = sep_dim == 0 ? left_child : this->sign_a;
        qpart.sign_y = sep_dim == 1 ? left_child : (sep_dim == 0 ? this->sign_a : this->sign_b);
        qpart.sign_z = sep_dim == 2 ? left_child : this->sign_b;
        return qpart;
    }
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol