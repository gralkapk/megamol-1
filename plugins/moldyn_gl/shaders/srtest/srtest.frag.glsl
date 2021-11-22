#version 450

flat in vec3 objPos;
flat in float rad;
flat in float sqrRad;
flat in vec4 pointColor;
flat in vec3 oc_pos;

#include "srtest_ubo.glsl"

layout(location = 0) out vec4 outColor;
layout(depth_greater) out float gl_FragDepth;

#include "lightdirectional.glsl"

#include "srtest_intersection.glsl"

#include "srtest_depth.glsl"

void main(void) {
    /*vec4 new_pos;
    vec3 normal;
    vec3 ray;
    float t;
    intersection(objPos, sqrRad, oc_pos, c, rad, new_pos, normal, ray, t);*/

    vec4 pos_ndc =
        vec4(2.0f * (gl_FragCoord.xy / viewAttr.zw) - 1.0f, (2.0f * gl_FragCoord.z) / (far - near) - 1.0f, 1.0f);
    vec4 pos_clip = MVPinv * pos_ndc;
    vec3 pos_obj = pos_clip.xyz / pos_clip.w;

    vec3 ray = normalize(pos_obj - camPos);

    float tf = dot(oc_pos, ray.xyz);
    vec3 tt = tf * ray.xyz - oc_pos;
    float delta = sqrRad - dot(tt, tt);
    if (delta < 0.0f)
        discard;

    float tb = sqrt(delta);
    float t = tf - tb;

    vec4 new_pos = vec4(camPos + t * ray.xyz, 1.0f);

    vec3 normal = (new_pos.xyz - objPos) / rad;

    outColor = vec4(LocalLighting(ray.xyz, normal, lightDir, pointColor.rgb), pointColor.a);
    /*outColor = vec4(factor, 0, 1);
    outColor = vec4(ve0, 1);*/
    //outColor = vec4(0.5f*(ray+1.0f), 1);

    gl_FragDepth = depth(t);
}
