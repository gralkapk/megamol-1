#version 450

uniform vec4 globalCol;
uniform float globalRad;

uniform bool useGlobalCol;
uniform bool useGlobalRad;

out VPoint {
    flat vec3  objPos;
    flat float rad;
    flat float sqrRad;
    flat vec4  pointColor;
    flat vec3  oc_pos;
}
v_pp;

#include "srtest_ubo.glsl"

#ifdef __SRTEST_VAO__
#include "srtest_vao.glsl"
#elif defined(__SRTEST_SSBO__)
#include "srtest_ssbo.glsl"
#endif

void main() {
    access_data(gl_VertexID, v_pp.objPos, v_pp.pointColor, v_pp.rad);

    v_pp.oc_pos = v_pp.objPos - camPos;
    v_pp.sqrRad = v_pp.rad * v_pp.rad;

    gl_Position = vec4(v_pp.objPos, 1.0f);
}