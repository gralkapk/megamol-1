/**
 * MegaMol
 * Copyright (c) 2023, MegaMol Dev Team
 * All rights reserved.
 */
#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <TimeTypes.h>

#include "CPUQuery.h"
#include "GLQuery.h"

namespace megamol::frontend_resources::performance {

enum class parent_type { CALL, USER_REGION, BUILTIN };

static constexpr const char* parent_type_string(parent_type parent) {
    switch (parent) {
    case parent_type::CALL:
        return "Call";
    case parent_type::USER_REGION:
        return "UserRegion";
    case parent_type::BUILTIN:
        return "BuiltIn";
    }
    return "unknown";
}

struct timer_region {
    time_point start = zero_time, end = zero_time;
    int64_t global_index = -1;
    frame_type frame;
    std::array<std::shared_ptr<AnyQuery>, 2> qids;
    bool finished = false;
};

struct basic_timer_config {
    std::string name = "unnamed";
    query_api api = query_api::CPU;
    user_index_type user_index = 0;
};

struct timer_config : public basic_timer_config {
    parent_type parent = parent_type::CALL;
    void* parent_pointer = nullptr;
    std::string comment;
};

struct timer_entry {
    // the user cannot fiddle with timers directly, this class needs to be asked
    handle_type handle = 0;
    query_api api = query_api::CPU;
    // local index inside one frame (if this region is touched multiple times per frame)
    uint32_t frame_index = 0;
    // user payload, used to track call indices, for example
    user_index_type user_index = 0;
    parent_type parent = parent_type::BUILTIN;
    time_point start, end, duration;
    int64_t global_index = 0;
    frame_type frame = 0;
};

struct frame_info {
    frame_type frame = 0;
    std::vector<timer_entry> entries;
};

class Itimer {
    friend class PerformanceManager;

public:
    Itimer(timer_config conf) : conf(std::move(conf)) {}
    virtual ~Itimer() = default;

    [[nodiscard]] const timer_config& get_conf() const {
        return conf;
    }
    [[nodiscard]] handle_type get_handle() const {
        return h;
    }
    [[nodiscard]] uint32_t get_region_count() const {
        return regions.size();
    }
    [[nodiscard]] bool is_finished(uint32_t index) const {
        return regions[index].finished;
    }
    [[nodiscard]] time_point get_start(uint32_t index) const {
        return regions[index].start;
    }
    [[nodiscard]] time_point get_end(uint32_t index) const {
        return regions[index].end;
    }
    [[nodiscard]] int64_t get_global_index(uint32_t index) const {
        return regions[index].global_index;
    }
    [[nodiscard]] frame_type get_frame(uint32_t index) const {
        return regions[index].frame;
    }
    [[nodiscard]] frame_type get_start_frame() const {
        return start_frame;
    }

    // hint: this is not for free, so don't call this all the time
    static std::string parent_name(const timer_config& conf);

    void set_comment(std::string comment) {
        conf.comment = comment;
    }

protected:
    // returns whether this is a new frame from what has been seen
    virtual bool start(frame_type frame);

    virtual void end();

    virtual void collect(frame_type frame) = 0;

    virtual void clear(frame_type frame) = 0;

    timer_config conf;
    //time_point last_start;
    std::vector<timer_region> regions;
    bool started = false;
    frame_type start_frame = std::numeric_limits<frame_type>::max();
    handle_type h = 0;
    // there can only be one PerformanceManager currently.
    inline static int64_t current_global_index = 0;
};

template<class Q>
class any_timer : public Itimer {
    static_assert(std::is_base_of_v<AnyQuery, Q>, "Q must inherit from AnyQuery");
public:
    any_timer(const timer_config& cfg) : Itimer(cfg) {}

    bool start(frame_type frame) override {
        const auto new_frame = Itimer::start(frame);

        if (first_frame_) {
            timer_region r{zero_time, zero_time, current_global_index++, frame, {nullptr, nullptr}, false};
            r.qids[0] = std::make_shared<Q>();
            r.qids[0]->Counter();
            regions.emplace_back(r);
            first_frame_ = false;
        }

        return new_frame;
    }

    void end() override {
        Itimer::end();

        regions.back().qids[1] = std::make_shared<Q>();
        regions.back().qids[1]->Counter();

        timer_region r{
            zero_time, zero_time, current_global_index++, regions.back().frame + 1, {nullptr, nullptr}, false};
        r.qids[0] = regions.back().qids[1];
        regions.emplace_back(r);
    }

    void collect(frame_type frame) override {
        for (auto& r : regions) {
            time_point start_time = zero_time, end_time = zero_time;
            bool start_ready = false, end_ready = false;
            if (r.start == zero_time) {
                if (r.qids[0])
                    start_time = r.qids[0]->GetNW();

                if (start_time != zero_time) {
                    r.start = start_time;
                    start_ready = true;
                }
            } else {
                start_ready = true;
            }
            if (r.end == zero_time) {
                if (r.qids[1])
                    end_time = r.qids[1]->GetNW();

                if (end_time != zero_time) {
                    r.end = end_time;
                    end_ready = true;
                }
            } else {
                end_ready = true;
            }
            r.finished = start_ready && end_ready;
        }
    }

    void clear(frame_type frame) override {
        regions.erase(
            std::remove_if(regions.begin(), regions.end(), [](auto const& r) { return r.finished; }), regions.end());
    }

private:
    bool first_frame_ = true;
};

using cpu_timer = any_timer<CPUQuery>;
using gl_timer = any_timer<GLQuery>;

} // namespace megamol::frontend_resources::performance
