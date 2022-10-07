/*
 * ProbeBillboardGlyphRenderTasks.h
 *
 * Copyright (C) 2019 by Universitaet Stuttgart (VISUS).
 * All rights reserved.
 */

#ifndef PROBE_BILLBOARD_GLYPH_RENDER_TASK_H_INCLUDED
#define PROBE_BILLBOARD_GLYPH_RENDER_TASK_H_INCLUDED

#include <typeindex>

#include "mesh_gl/BaseRenderTaskRenderer.h"

#include "probe/ProbeCollection.h"

#include <imgui.h>

namespace megamol {
namespace probe_gl {

class ProbeGlyphRenderer : public mesh_gl::BaseRenderTaskRenderer {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) {
        return "ProbeGlyphRenderer";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) {
        return "...";
    }

    ProbeGlyphRenderer();
    ~ProbeGlyphRenderer();

protected:
    /**
     * The get extents callback. The module should set the members of
     * 'call' to tell the caller the extents of its data (bounding boxes
     * and times).
     *
     * @param call The calling call.
     *
     * @return The return value of the function.
     */
    bool GetExtents(mmstd_gl::CallRender3DGL& call) override;

    void createMaterialCollection() override;
    void updateRenderTaskCollection(mmstd_gl::CallRender3DGL& call, bool force_update) override;

    bool saveGlyphsToImages(core::param::ParamSlot& slot);

private:
    core::CallerSlot m_transfer_function_Slot;

    core::CallerSlot m_probes_slot;

    core::CallerSlot m_event_slot;

    core::param::ParamSlot m_billboard_size_slot;

    core::param::ParamSlot m_rendering_mode_slot;

    core::param::ParamSlot m_use_interpolation_slot;

    core::param::ParamSlot m_show_canvas_slot;

    core::param::ParamSlot m_canvas_color_slot;

    core::param::ParamSlot m_save_glyphs_to_textures_slot;

    std::shared_ptr<glowl::Mesh> m_billboard_dummy_mesh;

    std::shared_ptr<glowl::Texture2D> m_transfer_function;

    std::array<float, 2> m_tf_range;

    ImGuiContext* m_imgui_context;

    /**
     *
     */

    struct PerFrameData {
        int use_interpolation;

        int show_canvas;

        int padding0;
        int padding1;

        std::array<float, 4> canvas_color;

        GLuint64 tf_texture_handle;
        float tf_min;
        float tf_max;
    };

    struct TexturedGlyphData {
        glm::vec4 position;
        GLuint64 texture_handle;
        float slice_idx;
        float scale;
    };

    struct GlyphVectorProbeData {
        glm::vec4 position;
        glm::vec4 probe_direction;
        float scale;

        int probe_id;
        int state;

        float sample_cnt;
        std::array<float, 4> samples[32];
    };

    struct GlyphScalarProbeData {
        glm::vec4 position;
        glm::vec4 probe_direction;
        float scale;

        float sample_cnt;
        float samples[32];

        int probe_id;
        int state;
    };

    struct GlyphScalarDistributionProbeData {
        glm::vec4 position;
        glm::vec4 probe_direction;
        float scale;

        int probe_id;
        int state;

        float sample_cnt;
        std::array<float, 4> samples[32];
    };

    struct GlyphClusterIDData {
        glm::vec4 position;
        glm::vec4 probe_direction;
        float scale;

        int probe_id;
        int state;

        float sample_cnt;

        int cluster_id;
        int total_cluster_cnt; // we have some space to spare per glyph so why not...
        int padding1;
        int padding2;
    };

    struct GlyphBaseProbeData {
        glm::vec4 position;
        glm::vec4 probe_direction;
        float scale;

        int probe_id;
        int state;

        float padding0;
    };

    bool m_show_glyphs;

    std::vector<std::pair<std::type_index, size_t>> m_type_index_map;

    std::vector<std::string> m_textured_glyph_identifiers;
    std::vector<std::string> m_vector_probe_glyph_identifiers;
    std::vector<std::string> m_scalar_probe_glyph_identifiers;
    std::vector<std::string> m_scalar_distribution_probe_glyph_identifiers;
    std::vector<std::string> m_clusterID_glyph_identifiers;
    std::vector<std::string> m_dial_glyph_background_identifiers;

    std::vector<TexturedGlyphData> m_textured_glyph_data;
    std::vector<GlyphVectorProbeData> m_vector_probe_glyph_data;
    std::vector<GlyphScalarProbeData> m_scalar_probe_glyph_data;
    std::vector<GlyphScalarDistributionProbeData> m_scalar_distribution_probe_glyph_data;
    std::vector<GlyphClusterIDData> m_clusterID_glyph_data;
    std::vector<GlyphBaseProbeData> m_dial_glyph_background_data;

    std::vector<glowl::DrawElementsCommand> m_textured_gylph_draw_commands;
    std::vector<glowl::DrawElementsCommand> m_vector_probe_gylph_draw_commands;
    std::vector<glowl::DrawElementsCommand> m_scalar_probe_gylph_draw_commands;
    std::vector<glowl::DrawElementsCommand> m_scalar_distribution_probe_gylph_draw_commands;
    std::vector<glowl::DrawElementsCommand> m_clusterID_gylph_draw_commands;
    std::vector<glowl::DrawElementsCommand> m_dial_glyph_background_draw_commands;

    struct HierachicalScalarAIOGlyphs {
        struct Glyph {
            std::string boxplot_submesh;
            std::string histrogram_submesh;
            std::string value_arcs_submesh;
            std::string scatterplot_submesh;

            glowl::DrawElementsCommand circular_grid_draw_command;
            glowl::DrawElementsCommand boxplot_draw_command;
            glowl::DrawElementsCommand histogram_draw_command;
            glowl::DrawElementsCommand value_arcs_draw_command;
            glowl::DrawElementsCommand scatterplot_draw_command;

            struct Data {
                glm::vec4 position;
                glm::vec4 probe_direction;

                float scale;

                int probe_id;
                int state;

                float sample_cnt;

                //  int sample_offset;
                //  int padding0;
                //  int padding1;
                //  int padding2;
            };
        };

        // Probably don't need to store sample values or depth on the GPU since everything is baked into vertices
        //struct PerFrameData {
        //    float sample_value;
        //    float sample_depth;
        //};

        static constexpr char circular_grid_shader_identifier[]  = "HSAIOG_CircGrid";
        static constexpr char boxplot_shader_identifier[]        = "HSAIOG_Boxplot";
        static constexpr char histogram_shader_identifier[]      = "HSAIOG_Histogram";
        static constexpr char value_arcs_shader_identifier[]     = "HSAIOG_ValueArcs";
        static constexpr char scatterplot_shader_identifier[]    = "HSAIOG_Scatterplot";

        static constexpr char circular_grid_submesh_identifier[] = "HSAIOG_CircGrid";

        std::vector<PerFrameData> per_frame_data;
        std::vector<Glyph::Data> per_glpyh_data;
        std::vector<Glyph> glyphs;

        static void createShaders(std::shared_ptr<mesh_gl::GPUMaterialCollection> material_collection);

        void addGlyph(probe::FloatDistributionProbe const& probe,
            std::shared_ptr<mesh_gl::GPUMeshCollection> mesh_collection);
    };

    bool addAllRenderTasks();

    void updateAllRenderTasks();

    template<typename ProbeType>
    TexturedGlyphData createTexturedGlyphData(
        ProbeType const& probe, int probe_id, GLuint64 texture_handle, float slice_idx, float scale);

    GlyphScalarProbeData createScalarProbeGlyphData(probe::FloatProbe const& probe, int probe_id, float scale);

    GlyphScalarDistributionProbeData createScalarDistributionProbeGlyphData(
        probe::FloatDistributionProbe const& probe, int probe_id, float scale);

    GlyphVectorProbeData createVectorProbeGlyphData(probe::Vec4Probe const& probe, int probe_id, float scale);

    GlyphBaseProbeData createBaseProbeGlyphData(probe::BaseProbe const& probe, int probe_id, float scale);

    GlyphClusterIDData createClusterIDGlyphData(probe::BaseProbe const& probe, int probe_id, float scale);
};

template<typename ProbeType>
inline ProbeGlyphRenderer::TexturedGlyphData ProbeGlyphRenderer::createTexturedGlyphData(
    ProbeType const& probe, int probe_id, GLuint64 texture_handle, float slice_idx, float scale) {
    TexturedGlyphData glyph_data;
    glyph_data.position = glm::vec4(probe.m_position[0] + probe.m_direction[0] * (probe.m_begin * 1.1f),
        probe.m_position[1] + probe.m_direction[1] * (probe.m_begin * 1.1f),
        probe.m_position[2] + probe.m_direction[2] * (probe.m_begin * 1.1f), 1.0f);
    glyph_data.texture_handle = texture_handle;
    glyph_data.slice_idx = slice_idx;
    glyph_data.scale = scale;

    // glyph_data.probe_id = probe_id;

    return glyph_data;
}

} // namespace probe_gl
} // namespace megamol


#endif // !PROBE_BILLBOARD_GLYPH_RENDER_TASK_H_INCLUDED