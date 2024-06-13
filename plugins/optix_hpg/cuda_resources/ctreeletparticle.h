#pragma once

#include "box.h"

namespace megamol {
namespace optix_hpg {
namespace device {
struct CPKDlet {
    box3f bounds;
    box3u32 fbounds;
    unsigned int begin, end;
    char offset[3];
    unsigned short prefix[3];
};

struct CPKDParticle {
    unsigned int dim : 2;
    unsigned int x : 10;
    unsigned int y : 10;
    unsigned int z : 10;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
