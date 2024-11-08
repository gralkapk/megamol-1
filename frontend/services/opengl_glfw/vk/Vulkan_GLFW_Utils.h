#pragma once

#include <iostream>

#include <VkBootstrap.h>

#include <GLFW/glfw3.h>

#include "mmcore/utility/log/Log.h"

namespace megamol::frontend {
// https://github.com/charles-lunarg/vk-bootstrap
inline bool init_vulkan(GLFWwindow* window, vkb::Instance& vkb_inst, VkSurfaceKHR& vkb_surface, vkb::Device& vkb_device) {
    vkb::InstanceBuilder builder;
    auto inst_ret =
        builder
            .set_app_name("MegaMol")
#ifdef DEBUG
            .request_validation_layers()
#endif
            .set_debug_callback(
                [](VkDebugUtilsMessageSeverityFlagBitsEXT messageSeverity, VkDebugUtilsMessageTypeFlagsEXT messageType,
                    const VkDebugUtilsMessengerCallbackDataEXT* pCallbackData, void* pUserData) -> VkBool32 {
                    auto severity = vkb::to_string_message_severity(messageSeverity);
                    auto type = vkb::to_string_message_type(messageType);
                    switch (messageSeverity) {
                    case VK_DEBUG_UTILS_MESSAGE_SEVERITY_INFO_BIT_EXT:
                        core::utility::log::Log::DefaultLog.WriteInfo("Vulkan: [%s] %s", type, pCallbackData->pMessage);
                        break;
                    case VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT:
                        core::utility::log::Log::DefaultLog.WriteWarn("Vulkan: [%s] %s", type, pCallbackData->pMessage);
                        break;
                    case VK_DEBUG_UTILS_MESSAGE_SEVERITY_VERBOSE_BIT_EXT:
                        core::utility::log::Log::DefaultLog.WriteWarn(
                            "Vulkan: [%s: VERBOSE] %s", type, pCallbackData->pMessage);
                        break;
                    case VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT:
                    default:
                        core::utility::log::Log::DefaultLog.WriteError(
                            "Vulkan: [%s] %s", type, pCallbackData->pMessage);
                    };
                    return VK_FALSE;
                })
            .build();

    if (!inst_ret) {
        core::utility::log::Log::DefaultLog.WriteError(
            "Failed to create Vulkan instance. Error: %s", inst_ret.error().message());
        return false;
    }
    vkb_inst = inst_ret.value();

    glfwCreateWindowSurface(vkb_inst, window, nullptr, &vkb_surface);

    vkb::PhysicalDeviceSelector selector{vkb_inst};
    auto phys_ret = selector.set_surface(vkb_surface)
                        .set_minimum_version(1, 1) // require a vulkan 1.1 capable device
                        .add_required_extension("VK_KHR_swapchain")
                        .select();
    if (!phys_ret) {
        core::utility::log::Log::DefaultLog.WriteError(
            "Failed to select Vulkan Physical Device. Error: %s", phys_ret.error().message());
        vkb::destroy_surface(vkb_inst, vkb_surface);
        vkb::destroy_instance(vkb_inst);
        return false;
    }

    vkb::DeviceBuilder device_builder{phys_ret.value()};
    // automatically propagate needed data from instance & physical device
    auto dev_ret = device_builder.build();
    if (!dev_ret) {
        core::utility::log::Log::DefaultLog.WriteError(
            "Failed to create Vulkan device. Error: %s", dev_ret.error().message());
        vkb::destroy_surface(vkb_inst, vkb_surface);
        vkb::destroy_instance(vkb_inst);
        return false;
    }
    vkb_device = dev_ret.value();

    return true;
}

bool create_swapchain(vkb::Device const& device, vkb::Swapchain& swapchain, int const width, int const height) {
    vkb::SwapchainBuilder swapchain_builder{device};
    auto swap_ret =
        swapchain_builder.set_desired_extent(width, height)
            .set_desired_format(vk::SurfaceFormatKHR(vk::Format::eR8G8B8A8Srgb, vk::ColorSpaceKHR::eSrgbNonlinear))
            .build();
    if (!swap_ret) {
        core::utility::log::Log::DefaultLog.WriteError(
            "Failed to create swapchain. Error: %s", swap_ret.error().message());
        return false;
    }
    swapchain = swap_ret.value();
    return true;
}
} // namespace megamol::frontend
