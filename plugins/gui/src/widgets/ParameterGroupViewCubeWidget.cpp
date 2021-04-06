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

megamol::gui::PickableCube::PickableCube(void)
    : image_up_arrow()
    , shader(nullptr) {
}


bool megamol::gui::PickableCube::Draw(unsigned int picking_id, int& inout_face_id, int& inout_orientation_id, int& out_hovered_face_id,
                                      int& out_hovered_orientation_id, const glm::vec4& cube_orientation, ManipVector& pending_manipulations) {

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
            "\n"
            "uniform mat4 rot_mx; \n"
            "uniform mat4 model_mx; \n"
            "uniform mat4 proj_mx; \n"
            "uniform int face_id; \n"
            "uniform int orientation_id; \n"
            "uniform int face_hover_id; \n"
            "uniform int orientation_hover_id; \n"
            "out vec2 tex_coord; \n"
            "flat out vec3 texture_color; \n"
            "flat out vec3 vertex_color; \n"
            "flat out int current_face_id; \n"
            "flat out int current_orientation_id; \n"
            "\n"
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
            "    const vec3 colors[6] = vec3[6](vec3(0.0, 0.0, 1.0), vec3(0.0, 1.0, 1.0), \n"
            "                                   vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 1.0), \n"
            "                                   vec3(0.0, 1.0, 0.0), vec3(1.0, 1.0, 0.0)); \n"
            "    // IDs and indices \n"
            "    float vertex_index            = float(gl_VertexID); \n"
            "    float mod_index               = vertex_index - (12.0 * floor(vertex_index/12.0)); \n"
            "    int current_face_index        = int(gl_VertexID / 12);        // in range [0-5] \n"
            "    current_face_id               = current_face_index;           // same indices as  in megamol::core::view::AbstractView3D::DefaultView \n"
            "    int current_orientation_index = int(floor(mod_index / 3.0));  // in range [0-3] \n"
            "    current_orientation_id        = current_orientation_index;    // same indices as  in megamol::core::view::AbstractView3D::DefaultOrientation \n"
            "    \n"
            "    // Vertex Color \n"
            "    vertex_color = colors[current_face_index] * 0.3; \n"
            "    if (face_id == current_face_id) { \n"
            "        vertex_color = colors[current_face_index] * (0.4 + (0.6 - 0.6 * (float(current_orientation_index) / 3.0))); \n"
            "    } \n"
            "    if ((face_hover_id == current_face_id) && (orientation_hover_id == current_orientation_id)) { \n"
            "        vertex_color = colors[current_face_index] + vec3(0.5, 0.5, 0.5); \n"
            "        vertex_color = clamp(vertex_color, vec3(0.0, 0.0, 0.0), vec3(1.0, 1.0, 1.0)); \n"
            "        vertex_color *= (0.5 + (0.5 - 0.5 * (float(current_orientation_index) / 3.0))); \n"
            "    } \n"
            "    \n"
            "    // Up Arrow Texture \n"
            "    texture_color = vertex_color * 0.5; \n"
            "    if ((mod_index == 0) || (mod_index == 4))      tex_coord = vec2(1.0, 0.0); \n"
            "    else if ((mod_index == 1) || (mod_index == 9)) tex_coord = vec2(0.0, 0.0); \n"
            "    else if ((mod_index == 3) || (mod_index == 7)) tex_coord = vec2(1.0, 1.0); \n"
            "    else if ((mod_index == 6) || (mod_index == 10)) tex_coord = vec2(0.0, 1.0); \n"
            "    else if ((mod_index == 2) || (mod_index == 5) || "
            "             (mod_index == 8) || (mod_index == 11)) tex_coord = vec2(0.5, 0.5); \n"
            "    \n"
            "    // Vertex Position \n"
            "    gl_Position = proj_mx * model_mx * rot_mx * vertices[gl_VertexID]; \n"
            "}";

        std::string fragment_src =
            "#version 130  \n"
            "#extension GL_ARB_explicit_attrib_location : require \n"
            "\n"
            "#define BIT_OFFSET_ID           7 \n"
            "#define BIT_OFFSET_FACE         2 \n"
            "#define BIT_OFFSET_ORIENTATION  0 \n"
            "\n"
            "#define FACE_FRONT                  0 \n"
            "#define FACE_BACK                   1 \n"
            "#define FACE_RIGHT                  2 \n"
            "#define FACE_LEFT                   3 \n"
            "#define FACE_TOP                    4 \n"
            "#define FACE_BOTTOM                 5 \n"
            "#define CORNER_TOP_LEFT_FRONT       6 \n"
            "#define CORNER_TOP_RIGHT_FRONT      7 \n"
            "#define CORNER_TOP_LEFT_BACK        8 \n"
            "#define CORNER_TOP_RIGHT_BACK       9 \n"
            "#define CORNER_BOTTOM_LEFT_FRONT    10 \n"
            "#define CORNER_BOTTOM_RIGHT_FRONT   11 \n"
            "#define CORNER_BOTTOM_LEFT_BACK     12 \n"
            "#define CORNER_BOTTOM_RIGHT_BACK    13 \n"
            "#define EDGE_TOP_FRONT              14 \n"
            "#define EDGE_TOP_LEFT               15 \n"
            "#define EDGE_TOP_RIGHT              16 \n"
            "#define EDGE_TOP_BACK               17 \n"
            "#define EDGE_BOTTOM_FRONT           18 \n"
            "#define EDGE_BOTTOM_LEFT            19 \n"
            "#define EDGE_BOTTOM_RIGHT           20 \n"
            "#define EDGE_BOTTOM_BACK            21 \n"
            "#define EDGE_FRONT_LEFT             22 \n"
            "#define EDGE_FRONT_RIGHT            23 \n"
            "#define EDGE_BACK_LEFT              24 \n"
            "#define EDGE_BACK_RIGHT             25 \n"
            "\n"
            "in vec2 tex_coord; \n"
            "flat in vec3 texture_color; \n"
            "flat in vec3 vertex_color; \n"
            "flat in int current_face_id; \n"
            "flat in int current_orientation_id; \n"
            "uniform sampler2D tex; \n"
            "uniform int picking_id; \n"
            "uniform int face_hover_id; \n"
            "layout(location = 0) out vec4 outFragColor; \n"
            "layout(location = 1) out vec2 outFragInfo; \n"
            "\n"
            "float supersample(in vec2 uv, float w, float alpha) { \n"
            "    return smoothstep(0.5 - w, 0.5 + w, alpha); \n"
            "} \n"
            "\n"
            "void main(void) { \n"
            "    outFragColor = vec4(vertex_color, 1.0); \n"
            "    \n"
            "    // Arrow Texture - Supersample 4 extra points \n"
            "    float alpha = texture(tex, tex_coord).a; \n"
            "    if (alpha > 0.0) { \n"
            "        float smootingEdge = fwidth(alpha); \n"
            "        float dscale = 0.354; // half of 1/sqrt2; you can play with this \n"
            "        vec2 duv = dscale * (dFdx(tex_coord) + dFdy(tex_coord)); \n"
            "        vec4 box = vec4(tex_coord-duv, tex_coord+duv); \n"
            "        float asum = supersample(box.xy, smootingEdge, alpha) \n"
            "                   + supersample(box.zw, smootingEdge, alpha) \n"
            "                   + supersample(box.xw, smootingEdge, alpha) \n"
            "                   + supersample(box.zy, smootingEdge, alpha); \n"
            "        alpha = (alpha + 0.5 * asum) / 3.0; \n"
            "        if (alpha > 0.0) { \n"
            "            outFragColor = mix(vec4(vertex_color, 1.0), vec4(texture_color, 1.0), alpha); \n"
            "        } \n"
            "    } \n"
            "    \n"
            "    int id = current_face_id; \n"
            "    const float de = 0.1; // must be in [0,1] \n"
            "    const float dc = 0.2;  // must be in [0,1] \n"
            "    bool highlight_corner = false; \n"
            "    bool highlight_edge   = false; \n"
            "    // ----- Corners ----- \n"
            "    if ((tex_coord.x > (1.0 - dc)) && (tex_coord.y > (1.0 - dc))) {  // ----- BOTTOM RIGHT ----- \n "
            "        if (current_face_id == FACE_FRONT)       id = CORNER_BOTTOM_RIGHT_FRONT; \n"
            "        else if (current_face_id == FACE_BACK)   id = CORNER_BOTTOM_LEFT_BACK; \n"
            "        else if (current_face_id == FACE_RIGHT)  id = CORNER_BOTTOM_RIGHT_BACK; \n"
            "        else if (current_face_id == FACE_LEFT)   id = CORNER_BOTTOM_LEFT_FRONT; \n"
            "        else if (current_face_id == FACE_TOP)    id = CORNER_TOP_RIGHT_FRONT; \n"
            "        else if (current_face_id == FACE_BOTTOM) id = CORNER_BOTTOM_RIGHT_BACK; \n"
            "        if (((current_face_id == FACE_FRONT)  && (face_hover_id == CORNER_BOTTOM_RIGHT_FRONT)) || \n"
            "           ((current_face_id == FACE_BACK)    && (face_hover_id == CORNER_BOTTOM_LEFT_BACK))   || \n"
            "           ((current_face_id == FACE_RIGHT)   && (face_hover_id == CORNER_BOTTOM_RIGHT_BACK))  || \n"
            "           ((current_face_id == FACE_LEFT)    && (face_hover_id == CORNER_BOTTOM_LEFT_FRONT))  || \n"
            "           ((current_face_id == FACE_TOP)     && (face_hover_id == CORNER_TOP_RIGHT_FRONT))    || \n"
            "           ((current_face_id == FACE_BOTTOM)  && (face_hover_id == CORNER_BOTTOM_RIGHT_BACK)))  highlight_corner = true; \n"
            "    } else if ((tex_coord.x < dc) && (tex_coord.y < dc)) {           // ----- TOP LEFT ----- \n"
            "        if (current_face_id == FACE_FRONT)       id = CORNER_TOP_LEFT_FRONT; \n"
            "        else if (current_face_id == FACE_BACK)   id = CORNER_TOP_RIGHT_BACK; \n"
            "        else if (current_face_id == FACE_RIGHT)  id = CORNER_TOP_RIGHT_FRONT; \n"
            "        else if (current_face_id == FACE_LEFT)   id = CORNER_TOP_LEFT_BACK; \n"
            "        else if (current_face_id == FACE_TOP)    id = CORNER_TOP_LEFT_BACK; \n"
            "        else if (current_face_id == FACE_BOTTOM) id = CORNER_BOTTOM_LEFT_FRONT; \n"
            "        if (((current_face_id == FACE_FRONT) && (face_hover_id == CORNER_TOP_LEFT_FRONT))  || \n"
            "           ((current_face_id == FACE_BACK)   && (face_hover_id == CORNER_TOP_RIGHT_BACK))  || \n"
            "           ((current_face_id == FACE_RIGHT)  && (face_hover_id == CORNER_TOP_RIGHT_FRONT)) || \n"
            "           ((current_face_id == FACE_LEFT)   && (face_hover_id == CORNER_TOP_LEFT_BACK))   || \n"
            "           ((current_face_id == FACE_TOP)    && (face_hover_id == CORNER_TOP_LEFT_BACK))   || \n"
            "           ((current_face_id == FACE_BOTTOM) && (face_hover_id == CORNER_BOTTOM_LEFT_FRONT))) highlight_corner = true; \n"
            "    } else if ((tex_coord.x < dc) && (tex_coord.y > (1.0 - dc))) {   // ----- BOTTOM LEFT ----- \n"
            "        if (current_face_id == FACE_FRONT)       id = CORNER_BOTTOM_LEFT_FRONT; \n"
            "        else if (current_face_id == FACE_BACK)   id = CORNER_BOTTOM_RIGHT_BACK; \n"
            "        else if (current_face_id == FACE_RIGHT)  id = CORNER_BOTTOM_RIGHT_FRONT; \n"
            "        else if (current_face_id == FACE_LEFT)   id = CORNER_BOTTOM_LEFT_BACK; \n"
            "        else if (current_face_id == FACE_TOP)    id = CORNER_TOP_LEFT_FRONT; \n"
            "        else if (current_face_id == FACE_BOTTOM) id = CORNER_BOTTOM_LEFT_BACK; \n"
            "        if (((current_face_id == FACE_FRONT) && (face_hover_id == CORNER_BOTTOM_LEFT_FRONT))  || \n"
            "           ((current_face_id == FACE_BACK)   && (face_hover_id == CORNER_BOTTOM_RIGHT_BACK))  || \n"
            "           ((current_face_id == FACE_RIGHT)  && (face_hover_id == CORNER_BOTTOM_RIGHT_FRONT)) || \n"
            "           ((current_face_id == FACE_LEFT)   && (face_hover_id == CORNER_BOTTOM_LEFT_BACK))   || \n"
            "           ((current_face_id == FACE_TOP)    && (face_hover_id == CORNER_TOP_LEFT_FRONT))     || \n"
            "           ((current_face_id == FACE_BOTTOM) && (face_hover_id == CORNER_BOTTOM_LEFT_BACK))) highlight_corner = true; \n"
            "    } else if ((tex_coord.x > (1.0 - dc)) && (tex_coord.y < dc)) {   // ----- TOP RIGHT -----\n"
            "        if (current_face_id == FACE_FRONT)       id = CORNER_TOP_RIGHT_FRONT; \n"
            "        else if (current_face_id == FACE_BACK)   id = CORNER_TOP_LEFT_BACK; \n"
            "        else if (current_face_id == FACE_RIGHT)  id = CORNER_TOP_RIGHT_BACK; \n"
            "        else if (current_face_id == FACE_LEFT)   id = CORNER_TOP_LEFT_FRONT; \n"
            "        else if (current_face_id == FACE_TOP)    id = CORNER_TOP_RIGHT_BACK; \n"
            "        else if (current_face_id == FACE_BOTTOM) id = CORNER_BOTTOM_RIGHT_FRONT; \n"
            "        if (((current_face_id == FACE_FRONT) && (face_hover_id == CORNER_TOP_RIGHT_FRONT)) || \n"
            "           ((current_face_id == FACE_BACK)   && (face_hover_id == CORNER_TOP_LEFT_BACK))   || \n"
            "           ((current_face_id == FACE_RIGHT)  && (face_hover_id == CORNER_TOP_RIGHT_BACK))  || \n"
            "           ((current_face_id == FACE_LEFT)   && (face_hover_id == CORNER_TOP_LEFT_FRONT))  || \n"
            "           ((current_face_id == FACE_TOP)    && (face_hover_id == CORNER_TOP_RIGHT_BACK))  || \n"
            "           ((current_face_id == FACE_BOTTOM) && (face_hover_id == CORNER_BOTTOM_RIGHT_FRONT))) highlight_corner = true; \n"
            "    // ----- Edges ----- \n"
            "    } else if (tex_coord.x > (1.0 - de)) {                           // ----- RIGHT ----- \n"
            "        if (current_face_id == FACE_FRONT)       id = EDGE_FRONT_RIGHT; \n"
            "        else if (current_face_id == FACE_BACK)   id = EDGE_BACK_LEFT; \n"
            "        else if (current_face_id == FACE_RIGHT)  id = EDGE_BACK_RIGHT; \n"
            "        else if (current_face_id == FACE_LEFT)   id = EDGE_FRONT_LEFT; \n"
            "        else if (current_face_id == FACE_TOP)    id = EDGE_TOP_RIGHT; \n"
            "        else if (current_face_id == FACE_BOTTOM) id = EDGE_BOTTOM_RIGHT; \n"
            "        if (((current_face_id == FACE_FRONT) && (face_hover_id == EDGE_FRONT_RIGHT)) || \n"
            "           ((current_face_id == FACE_BACK)   && (face_hover_id == EDGE_BACK_LEFT))   || \n"
            "           ((current_face_id == FACE_RIGHT)  && (face_hover_id == EDGE_BACK_RIGHT))  || \n"
            "           ((current_face_id == FACE_LEFT)   && (face_hover_id == EDGE_FRONT_LEFT))  || \n"
            "           ((current_face_id == FACE_TOP)    && (face_hover_id == EDGE_TOP_RIGHT))   || \n"
            "           ((current_face_id == FACE_BOTTOM) && (face_hover_id == EDGE_BOTTOM_RIGHT))) highlight_edge = true; \n"
            "    } else if (tex_coord.x < de) {                                   // ----- LEFT ----- \n"
            "        if (current_face_id == FACE_FRONT)       id = EDGE_FRONT_LEFT; \n"
            "        else if (current_face_id == FACE_BACK)   id = EDGE_BACK_RIGHT; \n"
            "        else if (current_face_id == FACE_RIGHT)  id = EDGE_FRONT_RIGHT; \n"
            "        else if (current_face_id == FACE_LEFT)   id = EDGE_BACK_LEFT; \n"
            "        else if (current_face_id == FACE_TOP)    id = EDGE_TOP_LEFT; \n"
            "        else if (current_face_id == FACE_BOTTOM) id = EDGE_BOTTOM_LEFT; \n"
            "        if (((current_face_id == FACE_FRONT) && (face_hover_id == EDGE_FRONT_LEFT))  || \n"
            "           ((current_face_id == FACE_BACK)   && (face_hover_id == EDGE_BACK_RIGHT))  || \n"
            "           ((current_face_id == FACE_RIGHT)  && (face_hover_id == EDGE_FRONT_RIGHT)) || \n"
            "           ((current_face_id == FACE_LEFT)   && (face_hover_id == EDGE_BACK_LEFT))   || \n"
            "           ((current_face_id == FACE_TOP)    && (face_hover_id == EDGE_TOP_LEFT))    || \n"
            "           ((current_face_id == FACE_BOTTOM) && (face_hover_id == EDGE_BOTTOM_LEFT))) highlight_edge = true; \n"
            "    } else if (tex_coord.y > (1.0 - de)) {                           // ----- BOTTOM ----- \n"
            "        if (current_face_id == FACE_FRONT)       id = EDGE_BOTTOM_FRONT; \n"
            "        else if (current_face_id == FACE_BACK)   id = EDGE_BOTTOM_BACK; \n"
            "        else if (current_face_id == FACE_RIGHT)  id = EDGE_BOTTOM_RIGHT; \n"
            "        else if (current_face_id == FACE_LEFT)   id = EDGE_BOTTOM_LEFT; \n"
            "        else if (current_face_id == FACE_TOP)    id = EDGE_TOP_FRONT; \n"
            "        else if (current_face_id == FACE_BOTTOM) id = EDGE_BOTTOM_BACK; \n"
            "        if (((current_face_id == FACE_FRONT) && (face_hover_id == EDGE_BOTTOM_FRONT)) || \n"
            "           ((current_face_id == FACE_BACK)   && (face_hover_id == EDGE_BOTTOM_BACK))  || \n"
            "           ((current_face_id == FACE_RIGHT)  && (face_hover_id == EDGE_BOTTOM_RIGHT)) || \n"
            "           ((current_face_id == FACE_LEFT)   && (face_hover_id == EDGE_BOTTOM_LEFT))  || \n"
            "           ((current_face_id == FACE_TOP)    && (face_hover_id == EDGE_TOP_FRONT))    || \n"
            "           ((current_face_id == FACE_BOTTOM) && (face_hover_id == EDGE_BOTTOM_BACK))) highlight_edge = true; \n"
            "    } else if (tex_coord.y < de) {                                   // ----- TOP ----- \n"
            "        if (current_face_id == FACE_FRONT)       id = EDGE_TOP_FRONT; \n"
            "        else if (current_face_id == FACE_BACK)   id = EDGE_TOP_BACK; \n"
            "        else if (current_face_id == FACE_RIGHT)  id = EDGE_TOP_RIGHT; \n"
            "        else if (current_face_id == FACE_LEFT)   id = EDGE_TOP_LEFT; \n"
            "        else if (current_face_id == FACE_TOP)    id = EDGE_TOP_BACK; \n"
            "        else if (current_face_id == FACE_BOTTOM) id = EDGE_BOTTOM_FRONT; \n"
            "        if (((current_face_id == FACE_FRONT) && (face_hover_id == EDGE_TOP_FRONT))  || \n"
            "           ((current_face_id == FACE_BACK)   && (face_hover_id == EDGE_TOP_BACK))   || \n"
            "           ((current_face_id == FACE_RIGHT)  && (face_hover_id == EDGE_TOP_RIGHT))  || \n"
            "           ((current_face_id == FACE_LEFT)   && (face_hover_id == EDGE_TOP_LEFT))   || \n"
            "           ((current_face_id == FACE_TOP)    && (face_hover_id == EDGE_TOP_BACK))   || \n"
            "           ((current_face_id == FACE_BOTTOM) && (face_hover_id == EDGE_BOTTOM_FRONT))) highlight_edge = true; \n"
            "    } \n"
            "    \n"
            "    if (highlight_corner) outFragColor = vec4(1.0, 1.0, 1.0, 1.0); \n"
            "    if (highlight_edge)   outFragColor = vec4(1.0, 1.0, 1.0, 1.0); \n"
            "    // Ecoded ID \n"
            "    int encoded_id = int((picking_id             << BIT_OFFSET_ID)   | \n"
            "                         (id                     << BIT_OFFSET_FACE) | \n"
            "                         (current_orientation_id << BIT_OFFSET_ORIENTATION)); \n"
            "    outFragInfo  = vec2(float(encoded_id), gl_FragCoord.z); \n"
            "} ";

        if (!megamol::core::view::RenderUtils::CreateShader(this->shader, vertex_src, fragment_src)) {
            return false;
        }
    }

    // Process pending manipulations ------------------------------------------
    const int BIT_OFFSET_ID          = 7;
    const int BIT_OFFSET_FACE        = 2;  // uses 5 bits
    const int BIT_OFFSET_ORIENTATION = 0;  // uses 2 bits

    out_hovered_face_id        = -1;
    out_hovered_orientation_id = -1;
    for (auto& manip : pending_manipulations) {
        // Check for right picking ID
        if (picking_id == (picking_id & (manip.obj_id >> BIT_OFFSET_ID))) {
            // Extract face and orientation ID
            int picked_face_id        = static_cast<int>((manip.obj_id >> BIT_OFFSET_FACE)        & 0b11111);
            int picked_orientation_id = static_cast<int>((manip.obj_id >> BIT_OFFSET_ORIENTATION) & 0b11);

            if (manip.type == InteractionType::SELECT) {
                inout_face_id = picked_face_id;
                inout_orientation_id = picked_orientation_id;
                selected = true;
                std::cout << sizeof(int) << std::endl;
            } else if (manip.type == InteractionType::HIGHLIGHT) {
                out_hovered_face_id = picked_face_id;
                out_hovered_orientation_id = picked_orientation_id;
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
    int y = viewport[3] - size - static_cast<int>(ImGui::GetFrameHeightWithSpacing());
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
    this->shader->setUniform("face_id", inout_face_id);
    this->shader->setUniform("orientation_id", inout_orientation_id);
    this->shader->setUniform("face_hover_id", out_hovered_face_id);
    this->shader->setUniform("orientation_hover_id", out_hovered_orientation_id);
    this->shader->setUniform("picking_id", static_cast<int>(picking_id));

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


bool megamol::gui::PickableTexture::Draw(unsigned int picking_id, int face_id, int& out_orientation_change, int& out_hovered_arrow_id,
                                         ManipVector& pending_manipulations) {

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
            "uniform int picking_id; \n"
            "out vec2 tex_coord; \n"
            "flat out int encoded_id; \n"
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
            "    encoded_id  = int(picking_id << int(gl_VertexID / 6)); \n"
            "    gl_Position = vertices[gl_VertexID]; \n"
            "    tex_coord   = texcoords[gl_VertexID]; \n"
            "}";

        std::string fragment_src = "#version 130 \n"
                                   "#extension GL_ARB_explicit_attrib_location : require \n"
                                   "in vec2 tex_coord; \n"
                                   "flat in int encoded_id; \n"
                                   "uniform int face_id; \n"       // XXX -> Only supported when in [0,5]!
                                   "uniform sampler2D tex; \n"
                                   "uniform vec3 color; \n"
                                   "layout(location = 0) out vec4 outFragColor; \n"
                                   "layout(location = 1) out vec2 outFragInfo; \n"
                                   "float supersample(in vec2 uv, float w, float alpha) { \n"
                                   "    return smoothstep(0.5 - w, 0.5 + w, alpha); \n"
                                   "} \n"
                                   "void main(void) { \n"
                                   "    outFragInfo  = vec2(float(encoded_id), gl_FragCoord.z); \n"
                                   "    // Same colors as in vertex shader of pickable cube \n"
                                   "    const vec3 colors[6] = vec3[6](vec3(0.0, 0.0, 1.0), vec3(0.0, 1.0, 1.0), \n"
                                   "                                   vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 1.0), \n"
                                   "                                   vec3(0.0, 1.0, 0.0), vec3(1.0, 1.0, 0.0)); \n"
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
                                   "    outFragColor = vec4(colors[face_id] * 0.75, alpha); \n"
                                   "} ";

        if (!megamol::core::view::RenderUtils::CreateShader(this->shader, vertex_src, fragment_src)) {
            return false;
        }
    }

    // Process pending manipulations ------------------------------------------
    out_hovered_arrow_id = 0;
    for (auto& manip : pending_manipulations) {
        int orientation_change = 0;
        if (picking_id == (manip.obj_id >> 0)) {
            orientation_change = -1;
        } else if (picking_id == (manip.obj_id >> 1)) {
            orientation_change = 1;
        }
        if (orientation_change != 0) {
            if (manip.type == InteractionType::SELECT) {
                out_orientation_change = orientation_change;
                selected = true;
            } else if (manip.type == InteractionType::HIGHLIGHT) {
                out_hovered_arrow_id = orientation_change;
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
    int y = viewport[3] - size - static_cast<int>(ImGui::GetFrameHeightWithSpacing());
    glViewport(x, y, size, size);

    this->shader->use();

    auto texture_id = this->image_rotation_arrow.GetTextureID();
    if (texture_id != 0) {
        glEnable(GL_TEXTURE_2D);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture_id);
        glUniform1i(this->shader->getUniformLocation("tex"), static_cast<GLint>(0));
    }

    this->shader->setUniform("face_id", face_id);
    this->shader->setUniform("picking_id", static_cast<int>(picking_id));

    // Arrow Color
    glm::vec3 color(0.6, 0.6, 0.6);
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
    this->name = "defaultView";
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

    return (
        param_cubeOrientation && param_defaultView && param_showCube && param_defaultOrientation && param_resetView);
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
            int hovered_face = -1;
            int hovered_orientation = -1;
            auto cube_picking_id = param_defaultView->UID(); // Using any parameter UID
            inout_picking_buffer->AddInteractionObject(
                cube_picking_id, this->cube_widget.GetInteractions(cube_picking_id));
            bool selected = this->cube_widget.Draw(cube_picking_id, default_view, default_orientation,
                                                   hovered_face, hovered_orientation, cube_orientation,
                                                   inout_picking_buffer->GetPendingManipulations());
            int hovered_arrow = 0;
            int orientation_change = 0;
            auto arrow_picking_id = param_defaultOrientation->UID(); // Using any other parameter UID
            inout_picking_buffer->AddInteractionObject(
                    arrow_picking_id, this->texture_widget.GetInteractions(arrow_picking_id));
            this->texture_widget.Draw(arrow_picking_id, default_view, orientation_change,
                                      hovered_arrow, inout_picking_buffer->GetPendingManipulations());
            if (orientation_change < 0) {
                default_orientation = ((default_orientation-1) < 0)?(3):(default_orientation-1);
            }
            else if (orientation_change > 0) {
                default_orientation = ((default_orientation+1) > 3)?(0):(default_orientation+1);
            }

            if (selected) {
                param_resetView->ForceSetValueDirty();
            }
            param_defaultOrientation->SetValue(default_orientation);
            param_defaultView->SetValue(default_view);

            // Tooltip
            std::string tooltip_text;
            if (hovered_face >= 0) {
                switch (hovered_face) {
                case (DefaultView_t::DEFAULTVIEW_FACE_FRONT):
                    tooltip_text += "[Front]";
                    break;
                case (DefaultView_t::DEFAULTVIEW_FACE_BACK):
                    tooltip_text += "[Back]";
                    break;
                case (DefaultView_t::DEFAULTVIEW_FACE_RIGHT):
                    tooltip_text += "[Right]";
                    break;
                case (DefaultView_t::DEFAULTVIEW_FACE_LEFT):
                    tooltip_text += "[Left]";
                    break;
                case (DefaultView_t::DEFAULTVIEW_FACE_TOP):
                    tooltip_text += "[Top]";
                    break;
                    case (DefaultView_t::DEFAULTVIEW_FACE_BOTTOM):
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
                case (DefaultOrientation_t::DEFAULTORIENTATION_TOP):
                    tooltip_text += "0 degree";
                    break;
                case (DefaultOrientation_t::DEFAULTORIENTATION_RIGHT):
                    tooltip_text += "90 degree";
                    break;
                case (DefaultOrientation_t::DEFAULTORIENTATION_BOTTOM):
                    tooltip_text += "180 degree";
                    break;
                case (DefaultOrientation_t::DEFAULTORIENTATION_LEFT):
                    tooltip_text += "270 degree";
                    break;
                default:
                    break;
                }
            }

            /* DEACTIVATED
            if (hovered_arrow < 0) {
                tooltip_text = "Rotate Right";
            } else if (hovered_arrow > 0) {
                tooltip_text = "Rotate Left";
            } */

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
