#version 450

uniform vec4 viewAttr;

uniform vec3 camIn;
uniform vec3 camUp;
uniform vec3 camRight;

uniform mat4 MVinv;
uniform mat4 MVtransp;
uniform mat4 MVP;
uniform mat4 MVPinv;
uniform mat4 MVPtransp;

uniform vec4 lightDir;

uniform vec4 inConsts1;
uniform sampler1D colTab;
uniform float lengthScale;
uniform float lengthFilter;
uniform uint flagsAvailable;

uniform vec4 clipDat;
uniform vec4 clipCol;

#define CONSTRAD inConsts1.x
#define MIN_COLV inConsts1.y
#define MAX_COLV inConsts1.z
#define COLTAB_SIZE inConsts1.w

in Point {
    // r, r^2, len
    flat vec3 rad;
    flat vec4 camPos;
    flat vec4 objPos;
    flat vec4 rotLightDir;

    flat vec3 rotMatT0;
    flat vec3 rotMatT1;
    flat vec3 rotMatT2;
    flat uint discardFrag;

    flat vec4 objColor;
}
pp;

layout(location = 0) out vec4 outColor;

#define DEPTH
#define CLIP

#include "lightdirectional.glsl"

void main(void) {
    vec4 coord;
    vec3 ray, tmp;
    const float maxLambda = 50000.0;

    if (pp.discardFrag > uint(0)) {
        discard;
    }

    // transform fragment coordinates from window coordinates to view coordinates.
    coord = gl_FragCoord * vec4(viewAttr.z, viewAttr.w, 2.0, 0.0) + vec4(-1.0, -1.0, -1.0, 1.0);

    // transform fragment coordinates from view coordinates to object coordinates.
    coord = MVPinv * coord;
    coord /= coord.w;
    coord -= pp.objPos; // ... and move


    // calc the viewing ray
    ray = (pp.rotMatT0 * coord.x) + (pp.rotMatT1 * coord.y) + (pp.rotMatT2 * coord.z);
    ray = normalize(ray - pp.camPos.xyz);

    vec4 cpos = pp.camPos + vec4(pp.rad.z, 0.0, 0.0, 0.0);

    // calculate the geometry-ray-intersection

    // arrow parameters
#define CYL_RAD (pp.rad.x * 0.5)
#define CYL_RAD_SQ (pp.rad.y * 0.25)
#define CYL_LEN (pp.rad.z * 2.0)
#define TIP_RAD pp.rad.x
#define TIP_LEN (pp.rad.z * 0.8)

    // super unoptimized cone code

    float coneF = TIP_RAD / TIP_LEN;
    coneF *= coneF;
    float coneA = coneF * ray.x * ray.x - ray.y * ray.y - ray.z * ray.z;
    float coneB = 2.0 * (coneF * ray.x * cpos.x - ray.y * cpos.y - ray.z * cpos.z);
    float coneC = coneF * cpos.x * cpos.x - cpos.y * cpos.y - cpos.z * cpos.z;

    float rDc = dot(ray.yz, cpos.yz);
    float rDr = dot(ray.yz, ray.yz);

    vec2 radicand =
        vec2((rDc * rDc) - (rDr * (dot(cpos.yz, cpos.yz) - CYL_RAD_SQ)), coneB * coneB - 4.0 * coneA * coneC);
    vec2 divisor = vec2(rDr, 2.0 * coneA);
    vec2 radix = sqrt(radicand);
    vec2 minusB = vec2(-rDc, -coneB);

    vec4 lambda = vec4((minusB.x - radix.x) / divisor.x, (minusB.y + radix.y) / divisor.y,
        (minusB.x + radix.x) / divisor.x, (minusB.y - radix.y) / divisor.y);

    bvec4 invalid = bvec4((divisor.x == 0.0) || (radicand.x < 0.0), (divisor.y == 0.0) || (radicand.y < 0.0),
        (divisor.x == 0.0) || (radicand.x < 0.0), (divisor.y == 0.0) || (radicand.y < 0.0));

    vec4 ix = cpos.xxxx + ray.xxxx * lambda;


    invalid.x = invalid.x || (ix.x < TIP_LEN) || (ix.x > CYL_LEN);
    invalid.y = invalid.y || (ix.y < 0.0) || (ix.y > TIP_LEN);
    invalid.z = invalid.z || !(((ix.z > TIP_LEN) || (ix.x > CYL_LEN)) && (ix.z < CYL_LEN));
    invalid.w = invalid.w || !((ix.w > 0.0) && (ix.w < TIP_LEN));

    if (invalid.x && invalid.y && invalid.z && invalid.w) {
#ifdef CLIP
        discard;
#endif // CLIP
    }

    vec3 intersection, color;
    vec3 normal = vec3(1.0, 0.0, 0.0);
    color = pp.objColor.rgb;

    if (!invalid.y) {
        invalid.xzw = bvec3(true, true, true);
        intersection = cpos.xyz + (ray * lambda.y);
        normal = normalize(vec3(-TIP_RAD / TIP_LEN, normalize(intersection.yz)));
        //        color = vec3(1.0, 0.0, 0.0);
    }
    if (!invalid.x) {
        invalid.zw = bvec2(true, true);
        intersection = cpos.xyz + (ray * lambda.x);
        normal = vec3(0.0, normalize(intersection.yz));
    }
    if (!invalid.z) {
        invalid.w = true;
        lambda.z = (CYL_LEN - cpos.x) / ray.x;
        intersection = cpos.xyz + (ray * lambda.z);
    }
    if (!invalid.w) {
        lambda.w = (TIP_LEN - cpos.x) / ray.x;
        intersection = cpos.xyz + (ray * lambda.w);
    }

    //color.r = 1.0 - intersection.x / CYL_LEN;
    //color.g = 0.0; // intersection.x / CYL_LEN;
    //color.b = intersection.x / CYL_LEN;

    // phong lighting with directional light
    //gl_FragColor = vec4(color.rgb, 1.0);
    outColor = vec4(LocalLighting(ray, normal, pp.rotLightDir.xyz, color), 1.0);
    outColor.rgb = (0.75 * outColor.rgb) + (0.25 * color);

    // calculate depth
#ifdef DEPTH
    intersection -= (ray * dot(ray, vec3(pp.rad.z, 0.0, 0.0)));
    tmp = intersection;
    intersection.x = dot(pp.rotMatT0, tmp.xyz);
    intersection.y = dot(pp.rotMatT1, tmp.xyz);
    intersection.z = dot(pp.rotMatT2, tmp.xyz);

    intersection += pp.objPos.xyz;

    vec4 Ding = vec4(intersection, 1.0);
    float depth = dot(MVPtransp[2], Ding);
    float depthW = dot(MVPtransp[3], Ding);
#ifndef CLIP
    if (invalid.x && invalid.y && invalid.z && invalid.w) {
        gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
        gl_FragDepth = 0.99999;
    } else {
#endif // CLIP
        gl_FragDepth = ((depth / depthW) + 1.0) * 0.5;
#ifndef CLIP
    }
#endif // CLIP

    //    gl_FragColor.rgb *= ;

#endif // DEPTH

#ifdef RETICLE
    coord = gl_FragCoord * vec4(viewAttr.z, viewAttr.w, 2.0, 0.0) + vec4(-1.0, -1.0, -1.0, 1.0);
    if (min(abs(coord.x - centerFragment.x), abs(coord.y - centerFragment.y)) < 0.002) {
        gl_FragColor.rgb += vec3(0.3, 0.3, 0.5);
    }
#endif // RETICLE
}
