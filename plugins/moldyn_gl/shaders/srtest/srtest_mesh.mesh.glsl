#version 450

#extension GL_NV_mesh_shader : enable

layout(local_size_x = 32) in;
layout(max_vertices = 32, max_primitives = 32, points) out;

#include "srtest_ubo.glsl"

uniform uint num_points;

uniform vec4 globalCol;
uniform float globalRad;

uniform bool useGlobalCol;
uniform bool useGlobalRad;

#include "srtest_ssbo.glsl"

#include "srtest_touchplane.glsl"

out Point {
    flat vec3 objPos;
    flat float rad;
    flat float sqrRad;
    flat vec4 pointColor;
    flat vec3 oc_pos;
    flat float c;
}
pp[];

void main() {
    uint g_idx = gl_GlobalInvocationID.x; // TODO Check size
    if (g_idx < num_points) {
        uint l_idx = gl_LocalInvocationID.x;

        access_data(g_idx, pp[l_idx].objPos, pp[l_idx].pointColor, pp[l_idx].rad);

        pp[l_idx].oc_pos = camPos - pp[l_idx].objPos;
        pp[l_idx].sqrRad = pp[l_idx].rad * pp[l_idx].rad;
        pp[l_idx].c = dot(pp[l_idx].oc_pos, pp[l_idx].oc_pos) - pp[l_idx].sqrRad;

        vec4 projPos;
        float l;
        touchplane(pp[l_idx].objPos, pp[l_idx].rad, projPos, l);

        gl_MeshVerticesNV[l_idx].gl_Position = projPos;
        gl_MeshVerticesNV[l_idx].gl_PointSize = l;

        gl_PrimitiveIndicesNV[l_idx] = l_idx;
    }
    gl_PrimitiveCountNV = min(num_points - gl_WorkGroupID.x * gl_WorkGroupSize.x, gl_WorkGroupSize.x);
}
