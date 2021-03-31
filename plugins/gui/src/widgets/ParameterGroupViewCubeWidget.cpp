/*
 * ParameterGroupViewCubeWidget.cpp
 *
 * Copyright (C) 2019 by Universitaet Stuttgart (VIS).
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "widgets/ParameterGroupViewCubeWidget.h"
#include "graph/ParameterGroups.h"


using namespace megamol;
using namespace megamol::core;
using namespace megamol::gui;


// *** Pickable Cube ******************************************************** //

megamol::gui::PickableCube::PickableCube(void) : image_up_arrow(),shader(nullptr) {}


bool megamol::gui::PickableCube::Draw(unsigned int id, int& inout_face_index, int& inout_orientation_index,
                                      int& out_hovered_face, int& out_hovered_orientation, const glm::vec4& cube_orientation, ManipVector& pending_manipulations) {

    assert(ImGui::GetCurrentContext() != nullptr);
    bool selected = false;

    // Create texture once
    if (!this->image_up_arrow.IsLoaded()) {
        this->image_up_arrow.LoadTextureFromFile(GUI_VIEWCUBE_UP_ARROW);
    }

    // Create shader once -----------------------------------------------------
    if (this->shader == nullptr) {
        // INFO: IDs of the six cube faces are encoded via bit shift by face index of given parameter id.
        std::string vertex_src =
            "#version 130 \n"
            "uniform int id; \n"
            "uniform mat4 rot_mx; \n"
            "uniform mat4 model_mx; \n"
            "uniform mat4 proj_mx; \n"
            "uniform int face_index; \n"
            "uniform int orientation_index; \n"
            "uniform int face_hover_index; \n"
            "uniform int orientation_hover_index; \n"
            "out vec2 tex_coord; \n"
            "flat out vec4 vertex_color; \n"
            "flat out int face_id; \n"
            "flat out int use_texture; \n"
            "void main() { \n"
            "    // Vertex indices must fit enum order in megamol::core::view::View3D_2::defaultview \n"
            "    const vec4 vertices[72] = vec4[72]( \n"
            "        // DEFAULTVIEW_FRONT = 0 \n"
            "        vec4(1.0, 1.0, 1.0, 1.0), vec4(-1.0, 1.0, 1.0, 1.0), vec4(0.0, 0.0, 1.0, 1.0), \n"
            "        vec4(1.0, -1.0, 1.0, 1.0), vec4(1.0, 1.0, 1.0, 1.0), vec4(0.0, 0.0, 1.0, 1.0), \n"
            "        vec4(-1.0, -1.0, 1.0, 1.0), vec4(1.0, -1.0, 1.0, 1.0), vec4(0.0, 0.0, 1.0, 1.0), \n"
            "        vec4(-1.0, 1.0, 1.0, 1.0), vec4(-1.0, -1.0, 1.0, 1.0), vec4(0.0, 0.0, 1.0, 1.0), \n"
            "        // DEFAULTVIEW_BACK = 1 \n"
            "        vec4(-1.0, 1.0, -1.0, 1.0), vec4(1.0, 1.0, -1.0, 1.0), vec4(0.0, 0.0, -1.0, 1.0), \n"
            "        vec4(-1.0, -1.0, -1.0, 1.0), vec4(-1.0, 1.0, -1.0, 1.0), vec4(0.0, 0.0, -1.0, 1.0), \n"
            "        vec4(1.0, -1.0, -1.0, 1.0), vec4(-1.0, -1.0, -1.0, 1.0), vec4(0.0, 0.0, -1.0, 1.0), \n"
            "        vec4(1.0, 1.0, -1.0, 1.0), vec4(1.0, -1.0, -1.0, 1.0), vec4(0.0, 0.0, -1.0, 1.0), \n"
            "        // DEFAULTVIEW_RIGHT = 2 \n"
            "        vec4(1.0, 1.0, -1.0, 1.0), vec4(1.0, 1.0, 1.0, 1.0), vec4(1.0, 0.0, 0.0, 1.0), \n"
            "        vec4(1.0, -1.0, -1.0, 1.0), vec4(1.0, 1.0, -1.0, 1.0), vec4(1.0, 0.0, 0.0, 1.0), \n"
            "        vec4(1.0, -1.0, 1.0, 1.0), vec4(1.0, -1.0, -1.0, 1.0), vec4(1.0, 0.0, 0.0, 1.0), \n"
            "        vec4(1.0, 1.0, 1.0, 1.0), vec4(1.0, -1.0, 1.0, 1.0), vec4(1.0, 0.0, 0.0, 1.0), \n"
            "        // DEFAULTVIEW_LEFT = 3 \n"
            "        vec4(-1.0, 1.0, 1.0, 1.0), vec4(-1.0, 1.0, -1.0, 1.0), vec4(-1.0, 0.0, 0.0, 1.0), \n"
            "        vec4(-1.0, -1.0, 1.0, 1.0), vec4(-1.0, 1.0, 1.0, 1.0), vec4(-1.0, 0.0, 0.0, 1.0), \n"
            "        vec4(-1.0, -1.0, -1.0, 1.0), vec4(-1.0, -1.0, 1.0, 1.0), vec4(-1.0, 0.0, 0.0, 1.0), \n"
            "        vec4(-1.0, 1.0, -1.0, 1.0), vec4(-1.0, -1.0, -1.0, 1.0), vec4(-1.0, 0.0, 0.0, 1.0), \n"
            "        // DEFAULTVIEW_TOP = 4 \n"
            "        vec4(1.0, 1.0, -1.0, 1.0), vec4(-1.0, 1.0, -1.0, 1.0), vec4(0.0, 1.0, 0.0, 1.0), \n"
            "        vec4(1.0, 1.0, 1.0, 1.0), vec4(1.0, 1.0, -1.0, 1.0), vec4(0.0, 1.0, 0.0, 1.0), \n"
            "        vec4(-1.0, 1.0, 1.0, 1.0), vec4(1.0, 1.0, 1.0, 1.0), vec4(0.0, 1.0, 0.0, 1.0), \n"
            "        vec4(-1.0, 1.0, -1.0, 1.0), vec4(-1.0, 1.0, 1.0, 1.0), vec4(0.0, 1.0, 0.0, 1.0), \n"
            "        // DEFAULTVIEW_BOTTOM = 5 \n"
            "        vec4(1.0, -1.0, 1.0, 1.0), vec4(-1.0, -1.0, 1.0, 1.0), vec4(0.0, -1.0, 0.0, 1.0), \n"
            "        vec4(1.0, -1.0, -1.0, 1.0), vec4(1.0, -1.0, 1.0, 1.0), vec4(0.0, -1.0, 0.0, 1.0), \n"
            "        vec4(-1.0, -1.0, -1.0, 1.0), vec4(1.0, -1.0, -1.0, 1.0), vec4(0.0, -1.0, 0.0, 1.0), \n"
            "        vec4(-1.0, -1.0, 1.0, 1.0), vec4(-1.0, -1.0, -1.0, 1.0), vec4(0.0, -1.0, 0.0, 1.0)); \n"
            "    \n"
            "    const vec4 colors[6] = vec4[6](vec4(0.0, 0.0, 1.0, 1.0), vec4(0.0, 1.0, 1.0, 1.0), \n"
            "                                   vec4(0.0, 1.0, 0.0, 1.0), vec4(1.0, 1.0, 0.0, 1.0), \n"
            "                                   vec4(1.0, 0.0, 0.0, 1.0), vec4(1.0, 0.0, 1.0, 1.0)); \n"
            "    // Calculate indices and IDs \n"
            "    float vertex_index = float(gl_VertexID); \n"
            "    float mod_index = vertex_index - (12.0 * floor(vertex_index/12.0)); \n"
            "    float mod_triangle_index = mod_index - (3.0 * floor(mod_index/3.0)); \n"
            "    int current_orientation_index = int(floor(mod_index / 3.0)); \n"
            "    int current_orientation_id = int(1 << current_orientation_index); // in range [0-3]\n"
            "    int current_face_index = int(gl_VertexID / 12);       // in range [0-5] \n"
            "    face_id = int((id << (current_face_index + 4)) | current_orientation_id);"
            "    \n"
            "    // Set colors depending on selected or hovered triangles \n"
            "    vertex_color = colors[current_face_index] * 0.25; \n"
            "    if (face_index == current_face_index) { \n"
            "        vertex_color = colors[current_face_index] * (0.5 + (0.5 - 0.5*(current_orientation_index/2.0))); \n"
            "    } \n"
            "    if ((face_hover_index == current_face_index) && \n"
            "            (orientation_hover_index == current_orientation_index)) { \n"
            "        vertex_color = colors[current_face_index] * (0.6 + (0.4 - 0.4*(current_orientation_index/2.0))); \n"
            "    } \n"
            "    vertex_color.a = 1.0; \n"
            "    \n"
            "    use_texture = ((face_index == current_face_index)?(1):(0)); \n"
            "    if (use_texture > 0) { \n"
            "        if ((mod_index == 0) || (mod_index == 4)) tex_coord = vec2(1.0, 0.0); \n"
            "        else if ((mod_index == 1) || (mod_index == 9)) tex_coord = vec2(0.0, 0.0); \n"
            "        else if ((mod_index == 3) || (mod_index == 7)) tex_coord = vec2(1.0, 1.0); \n"
            "        else if ((mod_index == 6) || (mod_index == 10)) tex_coord = vec2(0.0, 1.0); \n"
            "        else if ((mod_index == 2) || (mod_index == 5) || (mod_index == 8) || (mod_index == 11)) tex_coord = vec2(0.5, 0.5); \n"
            "    } \n"
            "    gl_Position = proj_mx * model_mx * rot_mx * vertices[gl_VertexID]; \n"
            "}";

        std::string fragment_src = "#version 130  \n"
                                   "#extension GL_ARB_explicit_attrib_location : require \n"
                                   "in vec2 tex_coord; \n"
                                   "flat in vec4 vertex_color; \n"
                                   "flat in int face_id; \n"
                                   "flat in int use_texture; \n"
                                   "uniform sampler2D tex; \n"
                                   "uniform vec3 color; \n"
                                   "layout(location = 0) out vec4 outFragColor; \n"
                                   "layout(location = 1) out vec2 outFragInfo; \n"
                                   "float supersample(in vec2 uv, float w, float alpha) { \n"
                                   "    return smoothstep(0.5 - w, 0.5 + w, alpha); \n"
                                   "} \n"
                                   "void main(void) { \n"
                                   "    outFragColor = vertex_color; \n"
                                   "    outFragInfo  = vec2(float(face_id), gl_FragCoord.z); \n"
                                   "    if (use_texture > 0) {"
                                   "        vec4 tex_color = texture(tex, tex_coord); \n"
                                   "        float alpha = tex_color.a; \n"
                                   "        // Supersample - 4 extra points \n"
                                   "        float smootingEdge = fwidth(alpha); \n"
                                   "        float dscale = 0.354; // half of 1/sqrt2; you can play with this \n"
                                   "        vec2 duv = dscale * (dFdx(tex_coord) + dFdy(tex_coord)); \n"
                                   "        vec4 box = vec4(tex_coord-duv, tex_coord+duv); \n"
                                   "        float asum = supersample(box.xy, smootingEdge, alpha) \n"
                                   "                   + supersample(box.zw, smootingEdge, alpha) \n"
                                   "                   + supersample(box.xw, smootingEdge, alpha) \n"
                                   "                   + supersample(box.zy, smootingEdge, alpha); \n"
                                   "        alpha = (alpha + 0.5 * asum) / 3.0; \n"
                                   "        if (alpha > 0.0) outFragColor = mix(vertex_color, vec4(color, 1.0), alpha); \n"
                                   "    } \n"
                                   "} ";

        if (!megamol::core::view::RenderUtils::CreateShader(this->shader, vertex_src, fragment_src)) {
            return false;
        }
    }

    // Process pending manipulations ------------------------------------------
    out_hovered_orientation = -1;
    out_hovered_face = -1;
    for (auto& manip : pending_manipulations) {

        // ID is shifted by at least 4 bits and at most 9 bits.
        // Leaving at least 23 bit for actual id (meaning max id can be ....?)
        /// Indices must fit enum order in megamol::core::view::View3D_2::defaultview
        int face_index = -1;
        if (id == (id & (manip.obj_id >> 4))) // DEFAULTVIEW_FRONT
            face_index = 0;
        else if (id == (id & (manip.obj_id >> 5))) // DEFAULTVIEW_BACK
            face_index = 1;
        else if (id == (id & (manip.obj_id >> 6))) // DEFAULTVIEW_RIGHT
            face_index = 2;
        else if (id == (id & (manip.obj_id >> 7))) // DEFAULTVIEW_LEFT
            face_index = 3;
        else if (id == (id & (manip.obj_id >> 8))) // DEFAULTVIEW_TOP
            face_index = 4;
        else if (id == (id & (manip.obj_id >> 9))) // DEFAULTVIEW_BOTTOM
            face_index = 5;

        // First 4 bit indicate currently hovered orientation
        // Orientation is given by triangle order in shader of pickable cube
        int orientation_index = -1;
        if ((1 << 0) & manip.obj_id)
            orientation_index = 0; // TOP
        else if ((1 << 1) & manip.obj_id)
            orientation_index = 1; // RIGHT
        else if ((1 << 2) & manip.obj_id)
            orientation_index = 2; // BOTTOM
        else if ((1 << 3) & manip.obj_id)
            orientation_index = 3; // LEFT

        if (face_index >= 0) {
            if (manip.type == InteractionType::SELECT) {
                inout_face_index = face_index;
                inout_orientation_index = orientation_index;
                selected = true;
            } else if (manip.type == InteractionType::HIGHLIGHT) {
                out_hovered_face = face_index;
                out_hovered_orientation = orientation_index;
            }
        }
    }

    // Draw -------------------------------------------------------------------

    // Create view/model and projection matrices
    const auto rotation = glm::inverse(
        glm::mat4_cast(glm::quat(cube_orientation.w, cube_orientation.x, cube_orientation.y, cube_orientation.z)));
    const float dist = 2.0f / std::tan(megamol::core::thecam::math::angle_deg2rad(30.0f) / 2.0f);
    glm::mat4 model(1.0f);
    model[3][2] = -dist;
    const auto proj = glm::perspective(megamol::core::thecam::math::angle_deg2rad(30.0f), 1.0f, 0.1f, 100.0f);

    // Set state
    const auto culling = glIsEnabled(GL_CULL_FACE);
    if (!culling) {
        glEnable(GL_CULL_FACE);
    }
    std::array<GLint, 4> viewport;
    glGetIntegerv(GL_VIEWPORT, viewport.data());
    int size = (100 * static_cast<int>(megamol::gui::gui_scaling.Get()));
    int x = viewport[2] - size;
    int y = viewport[3] - size - ImGui::GetFrameHeightWithSpacing();
    glViewport(x, y, size, size);

    this->shader->use();

    auto texture_id = this->image_up_arrow.GetTextureID();
    if (texture_id != 0) {
        glEnable(GL_TEXTURE_2D);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture_id);
        glUniform1i(this->shader->getUniformLocation("tex"), static_cast<GLint>(0));
    }

    this->shader->setUniform("rot_mx", rotation);
    this->shader->setUniform("model_mx", model);
    this->shader->setUniform("proj_mx", proj);
    this->shader->setUniform("face_index", inout_face_index);
    this->shader->setUniform("orientation_index", inout_orientation_index);
    this->shader->setUniform("face_hover_index", out_hovered_face);
    this->shader->setUniform("orientation_hover_index", out_hovered_orientation);
    this->shader->setUniform("id", static_cast<int>(id));

    // Arrow Color
    glm::vec3 color(0.6, 0.6, 0.7);
    this->shader->setUniform("color", color);

    glDrawArrays(GL_TRIANGLES, 0, 72);

    if (texture_id != 0) {
        glBindTexture(GL_TEXTURE_1D, 0);
        glDisable(GL_TEXTURE_2D);
    }

    glUseProgram(0);

    // Restore
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    if (!culling) {
        glDisable(GL_CULL_FACE);
    }

    return selected;
}


InteractVector megamol::gui::PickableCube::GetInteractions(unsigned int id) const {

    InteractVector interactions;
    interactions.emplace_back(Interaction({InteractionType::SELECT, id, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}));
    interactions.emplace_back(Interaction({InteractionType::HIGHLIGHT, id, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}));
    return interactions;
}



// *** Pickable Cube ******************************************************** //

megamol::gui::PickableTexture::PickableTexture(void) : image_rotation_arrow(), shader(nullptr) {}


bool megamol::gui::PickableTexture::Draw(unsigned int id, int& out_orientation_index_offset, int& out_hovered_arrow, ManipVector& pending_manipulations) {

    assert(ImGui::GetCurrentContext() != nullptr);
    bool selected = false;

    // Create texture once
    if (!this->image_rotation_arrow.IsLoaded()) {
        this->image_rotation_arrow.LoadTextureFromFile(GUI_VIEWCUBE_ROTATION_ARROW);
    }

    // Create shader once -----------------------------------------------------
    if (this->shader == nullptr) {
        std::string vertex_src =
            "#version 130 \n"
            "uniform int id; \n"
            "out vec2 tex_coord; \n"
            "flat out int arrow_id; \n"
            "void main() { \n"
            "    const vec4 vertices[12] = vec4[12]( \n"
            "        vec4(0.75, 0.75, -1.0, 1.0), vec4(1.0, 0.75, -1.0, 1.0), vec4(1.0, 1.0, -1.0, 1.0), \n"
            "        vec4(0.75, 0.75, -1.0, 1.0), vec4(1.0, 1.0, -1.0, 1.0), vec4(0.75, 1.0, -1.0, 1.0), \n"
            "        vec4(0.0, 0.75, -1.0, 1.0), vec4(0.25, 0.75, -1.0, 1.0), vec4(0.25, 1.0, -1.0, 1.0), \n"
            "        vec4(0.0, 0.75, -1.0, 1.0), vec4(0.25, 1.0, -1.0, 1.0), vec4(0.0, 1.0, -1.0, 1.0)); \n"
            "    const vec2 texcoords[12] = vec2[12]( \n"
            "        vec2(0.0, 1.0), vec2(1.0, 1.0), vec2(1.0, 0.0), \n"
            "        vec2(0.0, 1.0), vec2(1.0, 0.0), vec2(0.0, 0.0), \n"
            "        vec2(1.0, 1.0), vec2(0.0, 1.0), vec2(0.0, 0.0), \n"
            "        vec2(1.0, 1.0), vec2(0.0, 0.0), vec2(1.0, 0.0)); \n"
            "    arrow_id = int(id << int(gl_VertexID / 6)); \n"
            "    gl_Position = vertices[gl_VertexID]; \n"
            "    tex_coord = texcoords[gl_VertexID]; \n"
            "}";

        std::string fragment_src = "#version 130 \n"
                                   "#extension GL_ARB_explicit_attrib_location : require \n"
                                   "in vec2 tex_coord; \n"
                                   "flat in int arrow_id; \n"
                                   "uniform sampler2D tex; \n"
                                   "uniform vec3 color; \n"
                                   "layout(location = 0) out vec4 outFragColor; \n"
                                   "layout(location = 1) out vec2 outFragInfo; \n"
                                   "float supersample(in vec2 uv, float w, float alpha) { \n"
                                   "    return smoothstep(0.5 - w, 0.5 + w, alpha); \n"
                                   "} \n"
                                   "void main(void) { \n"
                                   "    vec4 tex_color = texture(tex, tex_coord); \n"
                                   "    float alpha = tex_color.a; \n"
                                   "    // Supersample - 4 extra points \n"
                                   "    float smootingEdge = fwidth(alpha); \n"
                                   "    float dscale = 0.354; // half of 1/sqrt2; you can play with this \n"
                                   "    vec2 duv = dscale * (dFdx(tex_coord) + dFdy(tex_coord)); \n"
                                   "    vec4 box = vec4(tex_coord-duv, tex_coord+duv); \n"
                                   "    float asum = supersample(box.xy, smootingEdge, alpha) \n"
                                   "               + supersample(box.zw, smootingEdge, alpha) \n"
                                   "               + supersample(box.xw, smootingEdge, alpha) \n"
                                   "               + supersample(box.zy, smootingEdge, alpha); \n"
                                   "    alpha = (alpha + 0.5 * asum) / 3.0; \n"
                                   "    if (alpha <= 0.0) discard; \n"
                                   "    outFragColor = vec4(color, alpha); \n"
                                   "    outFragInfo  = vec2(float(arrow_id), gl_FragCoord.z); \n"
                                   "} ";

        if (!megamol::core::view::RenderUtils::CreateShader(this->shader, vertex_src, fragment_src)) {
            return false;
        }
    }

    // Process pending manipulations ------------------------------------------
    out_hovered_arrow = 0;
    for (auto& manip : pending_manipulations) {
        int orientation_index_offset = 0;
        if (id == manip.obj_id) {
            orientation_index_offset = -1;
        }
        else if (id == (id & (manip.obj_id >> 1))) {
            orientation_index_offset = 1;
        }
        if (orientation_index_offset != 0) {
            if (manip.type == InteractionType::SELECT) {
                out_orientation_index_offset = orientation_index_offset;
                selected = true;
            } else if (manip.type == InteractionType::HIGHLIGHT) {
                out_hovered_arrow = orientation_index_offset;
            }
        }
    }

    // Draw -------------------------------------------------------------------

    // Set state
    const auto culling = glIsEnabled(GL_CULL_FACE);
    if (!culling) {
        glEnable(GL_CULL_FACE);
    }
    std::array<GLint, 4> viewport;
    glGetIntegerv(GL_VIEWPORT, viewport.data());
    int size = (2 * 100 * static_cast<int>(megamol::gui::gui_scaling.Get()));
    int x = viewport[2] - size;
    int y = viewport[3] - size - ImGui::GetFrameHeightWithSpacing();
    glViewport(x, y, size, size);

    this->shader->use();

    auto texture_id = this->image_rotation_arrow.GetTextureID();
    if (texture_id != 0) {
        glEnable(GL_TEXTURE_2D);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture_id);
        glUniform1i(this->shader->getUniformLocation("tex"), static_cast<GLint>(0));
    }

    this->shader->setUniform("id", static_cast<int>(id));

    // Arrow Color
    glm::vec3 color(0.6, 0.6, 0.7);
    this->shader->setUniform("color", color);

    glDrawArrays(GL_TRIANGLES, 0, 12);

    if (texture_id != 0) {
        glBindTexture(GL_TEXTURE_1D, 0);
        glDisable(GL_TEXTURE_2D);
    }

    glUseProgram(0);

    // Restore
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    if (!culling) {
        glDisable(GL_CULL_FACE);
    }

    return selected;
}


InteractVector megamol::gui::PickableTexture::GetInteractions(unsigned int id) const {

    InteractVector interactions;
    interactions.emplace_back(Interaction({InteractionType::SELECT, id, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}));
    interactions.emplace_back(Interaction({InteractionType::HIGHLIGHT, id, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}));
    return interactions;
}



// *** Parameter Group View Cube Widget ************************************ //

megamol::gui::ParameterGroupViewCubeWidget::ParameterGroupViewCubeWidget(void)
        : AbstractParameterGroupWidget(megamol::gui::GenerateUniqueID())
        , tooltip()
        , cube_widget()
        , texture_widget()
        , last_presentation(param::AbstractParamPresentation::Presentation::Basic) {

    this->InitPresentation(Param_t::GROUP_3D_CUBE);
    this->name = "view";
}


bool megamol::gui::ParameterGroupViewCubeWidget::Check(bool only_check, ParamPtrVector_t& params) {

    bool param_cubeOrientation = false;
    bool param_defaultView = false;
    bool param_defaultOrientation = false;
    bool param_resetView = false;
    bool param_showCube = false;
    for (auto& param_ptr : params) {
        if ((param_ptr->Name() == "cubeOrientation") && (param_ptr->Type() == Param_t::VECTOR4F)) {
            param_cubeOrientation = true;
        } else if ((param_ptr->Name() == "defaultView") && (param_ptr->Type() == Param_t::ENUM)) {
            param_defaultView = true;
        } else if ((param_ptr->Name() == "defaultOrientation") && (param_ptr->Type() == Param_t::ENUM)) {
            param_defaultOrientation = true;
        } else if ((param_ptr->Name() == "resetView") && (param_ptr->Type() == Param_t::BUTTON)) {
            param_resetView = true;
        } else if ((param_ptr->Name() == "showViewCube") && (param_ptr->Type() == Param_t::BOOL)) {
            param_showCube = true;
        }
    }

    return (param_cubeOrientation && param_defaultView && param_showCube && param_defaultOrientation && param_resetView);
}


bool megamol::gui::ParameterGroupViewCubeWidget::Draw(ParamPtrVector_t params, const std::string& in_module_fullname,
    const std::string& in_search, megamol::gui::Parameter::WidgetScope in_scope, PickingBuffer* inout_picking_buffer) {

    if (ImGui::GetCurrentContext() == nullptr) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "No ImGui context available. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        return false;
    }

    // Check required parameters ----------------------------------------------
    Parameter* param_cubeOrientation = nullptr;
    Parameter* param_defaultView = nullptr;
    Parameter* param_defaultOrientation = nullptr;
    Parameter* param_resetView = nullptr;
    Parameter* param_showCube = nullptr;
    /// Find specific parameters of group by name because parameter type can occure multiple times.
    for (auto& param_ptr : params) {
        if ((param_ptr->Name() == "cubeOrientation") && (param_ptr->Type() == Param_t::VECTOR4F)) {
            param_cubeOrientation = param_ptr;
        } else if ((param_ptr->Name() == "defaultView") && (param_ptr->Type() == Param_t::ENUM)) {
            param_defaultView = param_ptr;
        } else if ((param_ptr->Name() == "defaultOrientation") && (param_ptr->Type() == Param_t::ENUM)) {
            param_defaultOrientation = param_ptr;
        } else if ((param_ptr->Name() == "resetView") && (param_ptr->Type() == Param_t::BUTTON)) {
            param_resetView = param_ptr;
        } else if ((param_ptr->Name() == "showViewCube") && (param_ptr->Type() == Param_t::BOOL)) {
            param_showCube = param_ptr;
        }
    }
    if ((param_cubeOrientation == nullptr) || (param_defaultView == nullptr) || (param_defaultOrientation == nullptr) ||
        (param_showCube == nullptr)) {
        utility::log::Log::DefaultLog.WriteError("[GUI] Unable to find all required parameters by name "
                                                 "for '%s' group widget. [%s, %s, line %d]\n",
            this->name.c_str(), __FILE__, __FUNCTION__, __LINE__);
        return false;
    }

    // Parameter presentation -------------------------------------------------
    auto presentation = this->GetGUIPresentation();
    if (presentation != this->last_presentation) {
        param_showCube->SetValue((presentation == param::AbstractParamPresentation::Presentation::Group_3D_Cube));
        this->last_presentation = presentation;
    } else {
        if (std::get<bool>(param_showCube->GetValue())) {
            this->last_presentation = param::AbstractParamPresentation::Presentation::Group_3D_Cube;
            this->SetGUIPresentation(this->last_presentation);
        } else {
            this->last_presentation = param::AbstractParamPresentation::Presentation::Basic;
            this->SetGUIPresentation(this->last_presentation);
        }
    }

    if (presentation == param::AbstractParamPresentation::Presentation::Basic) {

        if (in_scope == Parameter::WidgetScope::LOCAL) {
            // LOCAL

            ParameterGroups::DrawGroupedParameters(
                this->name, params, in_module_fullname, in_search, in_scope, nullptr, nullptr, GUI_INVALID_ID);

            return true;

        } else if (in_scope == Parameter::WidgetScope::GLOBAL) {

            // no global implementation ...
            return true;
        }

    } else if (presentation == param::AbstractParamPresentation::Presentation::Group_3D_Cube) {

        if (in_scope == Parameter::WidgetScope::LOCAL) {
            // LOCAL

            ParameterGroups::DrawGroupedParameters(
                this->name, params, "", in_search, in_scope, nullptr, nullptr, GUI_INVALID_ID);

            return true;

        } else if (in_scope == Parameter::WidgetScope::GLOBAL) {
            // GLOBAL

            if (inout_picking_buffer == nullptr) {
                utility::log::Log::DefaultLog.WriteError(
                    "[GUI] Pointer to required picking buffer is nullptr. [%s, %s, line %d]\n", __FILE__, __FUNCTION__,
                    __LINE__);
                return false;
            }

            ImGui::PushID(this->uid);

            auto default_view = std::get<int>(param_defaultView->GetValue());
            auto default_orientation = std::get<int>(param_defaultOrientation->GetValue());
            auto cube_orientation = std::get<glm::vec4>(param_cubeOrientation->GetValue());
            int hovered_cube_face = -1;
            int hovered_orientation = -1;
            auto cube_picking_id = param_defaultView->UID();
            inout_picking_buffer->AddInteractionObject(cube_picking_id, this->cube_widget.GetInteractions(cube_picking_id));
            bool selected = this->cube_widget.Draw(cube_picking_id, default_view, default_orientation, hovered_cube_face,
                                                   hovered_orientation, cube_orientation, inout_picking_buffer->GetPendingManipulations());

            int hovered_arrow = 0;
            int orientation_index_offset = 0;
            auto rot_picking_id = param_defaultOrientation->UID();
            inout_picking_buffer->AddInteractionObject(rot_picking_id, this->texture_widget.GetInteractions(rot_picking_id));
            this->texture_widget.Draw(rot_picking_id, orientation_index_offset, hovered_arrow, inout_picking_buffer->GetPendingManipulations());

            default_orientation = (default_orientation + orientation_index_offset);
            default_orientation = (default_orientation < 0)?(3):(default_orientation % 4);
            if (selected) {
                param_resetView->ForceSetValueDirty();
            }
            param_defaultOrientation->SetValue(default_orientation);
            param_defaultView->SetValue(default_view);

            // Tooltip
            std::string tooltip_text;
            if (hovered_cube_face >= 0) {
                /// Indices must fit enum order in view::View3D_2::defaultview
                switch (hovered_cube_face) {
                case (0): // DEFAULTVIEW_FRONT
                    tooltip_text += "[Front]";
                    break;
                case (1): // DEFAULTVIEW_BACK
                    tooltip_text += "[Back]";
                    break;
                case (2): // DEFAULTVIEW_RIGHT
                    tooltip_text += "[Right]";
                    break;
                case (3): // DEFAULTVIEW_LEFT
                    tooltip_text += "[Left]";
                    break;
                case (4): // DEFAULTVIEW_TOP
                    tooltip_text += "[Top]";
                    break;
                case (5): // DEFAULTVIEW_BOTTOM
                    tooltip_text += "[Bottom]";
                    break;
                default:
                    break;
                }
            }
            // Order is given by triangle order in shader of pickable cube
            if (hovered_orientation >= 0) {
                tooltip_text += " ";
                switch (hovered_orientation) {
                case (0): // TOP
                    tooltip_text += "0 degree";
                    break;
                case (1): // RIGHT
                    tooltip_text += "90 degree";
                    break;
                case (2): // BOTTOM
                    tooltip_text += "180 degree";
                    break;
                case (3): // LEFT
                    tooltip_text += "270 degree";
                    break;
                default:
                    break;
                }
            }

            if (hovered_arrow < 0) {
                tooltip_text = "Rotate Right";
            }
            else if (hovered_arrow > 0) {
                tooltip_text = "Rotate Left";
            }

            if (!tooltip_text.empty()) {
                ImGui::BeginTooltip();
                ImGui::TextUnformatted(tooltip_text.c_str());
                ImGui::EndTooltip();
            }

            ImGui::PopID();

            return true;
        }
    }
    return false;
}
