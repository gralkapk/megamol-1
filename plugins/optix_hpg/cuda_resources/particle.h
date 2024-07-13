#pragma once

#include <array>

#include <glm/glm.hpp>

#include "optix/utils_host.h"

#include "qtreeletparticle.h"
#include "btreeletparticle.h"
#include "ctreeletparticle.h"

namespace megamol {
namespace optix_hpg {

typedef union {
    unsigned int ui;
    struct {
        unsigned char a;
        unsigned char b;
        unsigned char c;
        unsigned char d;
    } parts;
} byte_cast;

namespace device {
using color_t = glm::u8vec4;

struct Particle {
    //float x, y, z;
    glm::vec3 pos;
};

struct PKDParticle {
    glm::vec3 pos;
    float dim;
};

struct PKDlet {
    //! bounding box of all particles (including the radius)
    box3f bounds;
    //! begin/end range in the common particles array
    unsigned int begin, end;
};

struct FPKDLet {
    box3u32 bounds;
    unsigned int begin, end;
    char offset[3];
    unsigned short prefix[3];
};

//struct SPKDlet {
//    //! bounding box of all particles (including the radius)
//    box3f bounds;
//    //! begin/end range in the common particles array
//    size_t begin, end;
//    unsigned char sx, sy, sz;
//    glm::vec3 lower;
//};

static int const spkd_array_size = 3;
struct SPKDlet {
    box3f bounds;
    unsigned int begin, end;
    glm::vec3 lower;
    unsigned char sx[spkd_array_size], sy[spkd_array_size], sz[spkd_array_size];
};

struct QPKDParticle {
#if 0
    unsigned int dim : 2;
    unsigned int sx : 1;
    unsigned int x : 9;
    unsigned int sy : 1;
    unsigned int y : 9;
    unsigned int sz : 1;
    unsigned int z : 9;
#endif
#if 0
    unsigned int dim_x : 1;
    unsigned int x : 31;
    unsigned int dim_y : 1;
    unsigned int y : 31;
    unsigned int dim_z : 1;
    unsigned int z : 31;
#endif
    unsigned int dim : 2;
    unsigned int sign_x : 1;
    unsigned int sign_y : 1;
    unsigned int sign_z : 1;
    unsigned int x : 16;
    unsigned int y : 16;
    unsigned int z : 16;
};

struct FPKDParticle {
    int dim;
    glm::uvec3 pos;
};

struct SPKDParticle {
    unsigned char dim : 2;
    unsigned char sx_idx : 2;
    unsigned char sy_idx : 2;
    unsigned char sz_idx : 2;
    unsigned char x;
    unsigned char y;
    unsigned char z;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
