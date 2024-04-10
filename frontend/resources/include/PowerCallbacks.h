#pragma once

#include <filesystem>
#include <functional>
#include <string>

namespace megamol::frontend_resources {

static std::string PowerCallbacks_Req_Name = "PowerCallbacks";

struct PowerCallbacks {
    std::function<unsigned long()> signal_high;
    std::function<unsigned long()> signal_low;
    std::function<void()> signal_frame;
    std::function<void(std::string const&, std::string const&)> add_meta_key_value;
    std::function<std::filesystem::path()> get_output_path;
};

} // namespace megamol::frontend_resources
