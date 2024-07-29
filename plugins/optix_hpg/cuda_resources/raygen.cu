#include "camera.h"
#include "raygen.h"

#include "optix/random.h"
#include "optix/random_owl.h"
#include "optix/utils_device.h"

namespace megamol {
namespace optix_hpg {
namespace device {
// OptiX SDK
// Path tracer example

//
// Copyright (c) 2019, NVIDIA CORPORATION. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of NVIDIA CORPORATION nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

//
// Modified 2021 MegaMol Dev Team
//

// code partially from: https://github.com/UniStuttgart-VISUS/rtxpkd_ldav2020
// ======================================================================== //
// Copyright 2018-2019 Ingo Wald                                            //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

// ======================================================================== //
// Modified 2019-2020 VISUS - University of Stuttgart                       //
// ======================================================================== //


//#define RANDVEC3F glm::vec3(rnd(42), rnd(42), rnd(42))
#define RANDVEC3F glm::vec3(rnd(), rnd(), rnd())

inline __device__ glm::vec3 random_in_unit_sphere(owl::common::LCG<16>& rnd) {
    glm::vec3 p;
    do {
        p = 2.0f * RANDVEC3F - glm::vec3(1, 1, 1);
    } while (glm::dot(p, p) >= 1.0f);
    return p;
}

inline __device__ glm::vec4 traceRay(const RayGenData& self, Ray& ray, owl::common::LCG<16>& rnd /*, Random& rnd*/,
    PerRayData& prd, glm::vec4& bg, int maxBounces) {

    unsigned int p0 = 0;
    unsigned int p1 = 0;
    packPointer(&prd, p0, p1);

    glm::vec3 col(1.f);

    for (int depth = 0; true; ++depth) {
        prd.particleID = -1;

        optixTrace(self.world, (const float3&) ray.origin, (const float3&) ray.direction, ray.tmin, ray.tmax, 0,
            (OptixVisibilityMask) -1,
            /*rayFlags     */ OPTIX_RAY_FLAG_DISABLE_ANYHIT,
            /*SBToffset    */ 0,
            /*SBTstride    */ 1,
            /*missSBTIndex */ 0, p0, p1);
        if (prd.particleID == -1) {
            return glm::vec4(col * glm::vec3(0.8f), 1.0f);
        }

        glm::vec3 N = (ray.origin + prd.t * ray.direction) - prd.pos;
        if (glm::dot(N, ray.direction) > 0.f)
            N = -N;
        N = glm::normalize(N);

        if (maxBounces == 0) {
            return glm::vec4(prd.albedo * (.2f + .6f * fabsf(glm::dot(N, ray.direction))), 1.0f);
        }

        col *= prd.albedo;

        if (depth >= maxBounces)
            return glm::vec4(0.1f, 0.1f, 0.1f, 1.0f);

        auto scattered_origin = ray.origin + prd.t * ray.direction;
        auto scattered_direction = N + random_in_unit_sphere(rnd);
        ray = Ray(/* origin   : */ scattered_origin,
            /* direction: */ glm::normalize(scattered_direction),
            /* tmin     : */ 1e-3f,
            /* tmax     : */ 1e+8f);
    }
}

MM_OPTIX_RAYGEN_KERNEL(raygen_program)() {
    // printf("RAYGEN1\n");
    const RayGenData& self = getProgramData<RayGenData>();
    auto const index = optixGetLaunchIndex();
    auto const dim = optixGetLaunchDimensions();
    glm::ivec2 pixelID = glm::ivec2(index.x, index.y);

    if (pixelID.x >= self.fbSize.x)
        return;
    if (pixelID.y >= self.fbSize.y)
        return;
    const int pixelIdx = pixelID.x + self.fbSize.x * pixelID.y;
    const int pixel_index = pixelID.y * dim.x + pixelID.x;

    const FrameState* fs = &self.frameStateBuffer[0];

    /*auto frame_idx = self.colorBufferPtr[pixelIdx].w;
    if (fs->changed) {
        frame_idx = 0.0f;
        self.colorBufferPtr[pixelIdx].w = 0.0f;
    }*/
    // auto const old_col = self.colorBufferPtr[pixelIdx];

    /*float4 old_col;
    surf2Dread(&old_col, self.col_surf, pixelID.x * sizeof(float4), pixelID.y, cudaBoundaryModeZero);*/

    //unsigned int seed = tea<16>(pixelID.y * self.fbSize.x + pixelID.x, fs->frameIdx);

    owl::common::LCG<16> rnd_owl(pixel_index, fs->frameIdx);


    glm::vec4 col(0.f);
    glm::vec4 bg = fs->background;

    // printf("RAYGEN FS %f\n", fs->near);

    //auto i = fs->samplesPerPixel;

    float depth = FLT_MAX;

    PerRayData prd;
    for (int s = 0; s < fs->samplesPerPixel; ++s) {
        prd.countDepth = true;
        prd.ray_depth = FLT_MAX;
        float u = -fs->rw + (fs->rw + fs->rw) * (float(pixelID.x) + rnd_owl()) / self.fbSize.x;
        float v = -(fs->th + (-fs->th - fs->th) * (float(pixelID.y) + rnd_owl()) / self.fbSize.y);
        /*float u = -fs->rw + (fs->rw + fs->rw) * float(pixelID.x) / self.fbSize.x;
        float v = -(fs->th + (-fs->th - fs->th) * float(pixelID.y) / self.fbSize.y);*/
        auto ray = generateRay(*fs, u, v);
        col += traceRay(self, ray, rnd_owl, prd, bg, fs->maxBounces);
        depth = fminf(depth, prd.ray_depth);
    }

    col /= (float) fs->samplesPerPixel;
    // col.w = frame_idx + 1;
    //++col.w;

    //if (fs->frameIdx > 0) {
    //    const float a = 1.0f / static_cast<float>(fs->frameIdx + 1);
    //    col = lerp(glm::vec4(static_cast<float>(old_col.x), static_cast<float>(old_col.y),
    //                   static_cast<float>(old_col.z), static_cast<float>(old_col.w)),
    //        col, a);
    //    // col.w = frame_idx + 1;
    //}

    if (depth < FLT_MAX) {
        depth = (fs->depth_params.z / depth) - (fs->depth_params.x);
        depth = 0.5f * (depth + 1.0f);
    } else {
        depth = 1.f;
        col = bg;
    }
    //surf2Dwrite(depth, self.depth_surf, pixelID.x * sizeof(float), pixelID.y, cudaBoundaryModeZero);


    if (fs->frameIdx > 0 && fs->accumulate)
        col += self.accumBuffer[pixelIdx];
    self.accumBuffer[pixelIdx] = col;

    if (fs->accumulate)
        col /= float(fs->frameIdx + 1);
    surf2Dwrite(make_float4(col.r, col.g, col.b, col.a), self.col_surf, pixelID.x * sizeof(float4), pixelID.y,
        cudaBoundaryModeZero);
    /*surf2Dwrite(make_float4(1, 1, 1, 1), self.col_surf, pixelID.x * sizeof(float4), pixelID.y,
        cudaBoundaryModeZero);*/


    //surf2Dwrite(1, self.depth_surf, pixelID.x * sizeof(float), pixelID.y, cudaBoundaryModeZero);
}
} // namespace device
} // namespace optix_hpg
} // namespace megamol
