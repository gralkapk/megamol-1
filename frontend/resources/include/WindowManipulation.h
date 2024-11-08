/**
 * MegaMol
 * Copyright (c) 2020, MegaMol Dev Team
 * All rights reserved.
 */

#pragma once

#include "FrameStatsCallbacks.h"
#include "GL_STUB.h"

namespace megamol::frontend_resources {

static std::string WindowManipulation_Req_Name = "WindowManipulation";

struct WindowManipulation {
    void set_window_title(const char* title) const;
    void set_framebuffer_size(const unsigned int width, const unsigned int height) const;
    void set_window_position(const unsigned int width, const unsigned int height) const;
    void set_swap_interval(const unsigned int wait_frames) const; // DANGER: assumes there is a GL context active
    std::function<void(const int)> set_mouse_cursor;

    void swap_buffers() const;

    enum class Fullscreen { Maximize, Restore };
    void set_fullscreen(const Fullscreen action) const;

    void* window_ptr = nullptr;
};

} // namespace megamol::frontend_resources
