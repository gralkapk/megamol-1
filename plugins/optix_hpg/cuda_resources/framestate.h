#pragma once

#include <glm/glm.hpp>

#include "optix/Ray.h"

#undef near

namespace megamol {
namespace optix_hpg {
namespace device {
struct FrameState {
    glm::vec3 camera_up;
    glm::vec3 camera_right;
    glm::vec3 camera_front;
    glm::vec3 camera_center;

    float th;
    float rw;
    float near;

    int samplesPerPixel{1};
    int maxBounces{0};
    bool accumulate;

    glm::vec4 background;

    int frameIdx;

    glm::vec3 depth_params;

    float intensity;
};

struct FrameState_peel {
    glm::vec3 camera_up;
    glm::vec3 camera_right;
    glm::vec3 camera_front;
    glm::vec3 camera_center;

    float th;
    float rw;
    float near;

    int samplesPerPixel{1};
    int maxBounces{0};
    bool accumulate;

    glm::vec4 background;

    int frameIdx;

    glm::vec3 depth_params;

    float intensity;

    Ray* rayBuffer;
    glm::ivec2* treeletRequests;
};
} // namespace device
} // namespace optix_hpg
} // namespace megamol
