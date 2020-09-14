#pragma once

#include <unordered_map>
#include <functional>

namespace megamol {
namespace core {

class MegaMolGraph_Convenience {
public:
    // parameter groups hold requests for new parameter values
    // values are applied to the graph upon request at once
    struct ParameterGroup {
        std::string name;
        std::unordered_map<std::string, std::string> parameter_values;
        void* graph; // graph to apply queued values to

        bool QueueParameterValue(const std::string& id, const std::string& value);
        bool ApplyQueuedParameterValues();
    };

    MegaMolGraph_Convenience(void* graph_ptr = nullptr);

    ParameterGroup& CreateParameterGroup(const std::string& group_name);
    ParameterGroup* FindParameterGroup(const std::string& group_name);
    std::vector<std::reference_wrapper<ParameterGroup>> ListParameterGroups();

    bool CreateChainCall(const std::string callName, const std::string from_slot_name, const std::string to_slot_name);

private:
    void* m_graph_ptr = nullptr;

    std::unordered_map<std::string, ParameterGroup> m_parameter_groups;
};


} /* namespace core */
} // namespace megamol

