/*
 * ProbeShellElementsRenderTasks.h
 *
 * Copyright (C) 2021 by Universitaet Stuttgart (VISUS).
 * All rights reserved.
 */

#ifndef PROBE_SHELL_ELEMENTS_RENDER_TASK_H_INCLUDED
#define PROBE_SHELL_ELEMENTS_RENDER_TASK_H_INCLUDED

#include "mesh_gl/BaseMeshRenderer.h"

namespace megamol {
namespace probe_gl {

class ProbeShellElementsRenderTasks : public mesh_gl::BaseMeshRenderer {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) {
        return "ProbeShellElementsRenderer";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) {
        return "Mesh viewer for the elements of shells used for probing.";
    }

    ProbeShellElementsRenderTasks();
    ~ProbeShellElementsRenderTasks();

protected:
    void createMaterialCollection() override;
    void createRenderTaskCollection() override;
    void updateRenderTaskCollection(mmstd_gl::CallRender3DGL& call, bool force_update) override;

private:
    bool m_show_elements;

    std::vector<std::vector<std::string>> m_rt_identifiers;
    std::vector<std::shared_ptr<glowl::Mesh>> m_batch_meshes;
    std::vector<std::vector<glowl::DrawElementsCommand>> m_draw_commands;

    struct PerObjectData {
        std::array<float, 16> transform;
        std::array<float, 4> color;
        int cluster_id;
        int total_cluster_cnt;
        int padding1;
        int padding2;
    };
    std::vector<std::vector<PerObjectData>> m_per_object_data;

    glm::vec4 m_hull_color;

    core::CallerSlot m_probes_slot;

    core::CallerSlot m_event_slot;

    /** Slot for setting different rendering mode in hull shader */
    core::param::ParamSlot m_shading_mode_slot;

    /** Slot for setting the color of the hull mesh */
    core::param::ParamSlot m_hull_color_slot;
};

} // namespace probe_gl
} // namespace megamol


#endif // !PROBE_SHELL_ELEMENTS_RENDER_TASK_H_INCLUDED