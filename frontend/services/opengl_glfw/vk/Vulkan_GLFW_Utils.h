#pragma once

#include <iostream>

#include <VkBootstrap.h>

#include <GLFW/glfw3.h>

namespace megamol::frontend {
// https://github.com/charles-lunarg/vk-bootstrap
inline bool init_vulkan(GLFWwindow* window) {
    vkb::InstanceBuilder builder;
    auto inst_ret = builder.set_app_name("MegaMol").request_validation_layers().use_default_debug_messenger().build();

    if (!inst_ret) {
        std::cerr << "Failed to create Vulkan instance. Error: " << inst_ret.error().message() << "\n";
        return false;
    }
    vkb::Instance vkb_inst = inst_ret.value();

    VkSurfaceKHR surface;
    glfwCreateWindowSurface(vkb_inst, window, nullptr, &surface);

    vkb::PhysicalDeviceSelector selector{vkb_inst};
    auto phys_ret = selector.set_surface(surface)
                        .set_minimum_version(1, 1) // require a vulkan 1.1 capable device
                        .require_dedicated_transfer_queue()
                        .select();
    if (!phys_ret) {
        std::cerr << "Failed to select Vulkan Physical Device. Error: " << phys_ret.error().message() << "\n";
        return false;
    }

    vkb::DeviceBuilder device_builder{phys_ret.value()};
    // automatically propagate needed data from instance & physical device
    auto dev_ret = device_builder.build();
    if (!dev_ret) {
        std::cerr << "Failed to create Vulkan device. Error: " << dev_ret.error().message() << "\n";
        return false;
    }
    vkb::Device vkb_device = dev_ret.value();

    // Get the VkDevice handle used in the rest of a vulkan application
    VkDevice device = vkb_device.device;

    // Get the graphics queue with a helper function
    auto graphics_queue_ret = vkb_device.get_queue(vkb::QueueType::graphics);
    if (!graphics_queue_ret) {
        std::cerr << "Failed to get graphics queue. Error: " << graphics_queue_ret.error().message() << "\n";
        return false;
    }
    VkQueue graphics_queue = graphics_queue_ret.value();

    return true;
}
} // namespace megamol::frontend
