#pragma once

#include "glm/glm.hpp"

#define MM_OPTIX_RAYGEN_ANNOTATION __raygen__

#define MM_OPTIX_INTERSECTION_ANNOTATION __intersection__

#define MM_OPTIX_ANYHIT_ANNOTATION __anyhit__

#define MM_OPTIX_CLOSESTHIT_ANNOTATION __closesthit__

#define MM_OPTIX_MISS_ANNOTATION __miss__

#define MM_OPTIX_DIRECT_CALLABLE_ANNOTATION __direct_callable__

#define MM_OPTIX_CONTINUATION_CALLABLE_ANNOTATION __continuation_callable__

#define MM_OPTIX_EXCEPTION_ANNOTATION __exception__

#define MM_OPTIX_BOUNDS_ANNOTATION __boundsKernel__


#define MM_OPTIX_RAYGEN_ANNOTATION_STRING "__raygen__"

#define MM_OPTIX_INTERSECTION_ANNOTATION_STRING "__intersection__"

#define MM_OPTIX_ANYHIT_ANNOTATION_STRING "__anyhit__"

#define MM_OPTIX_CLOSESTHIT_ANNOTATION_STRING "__closesthit__"

#define MM_OPTIX_MISS_ANNOTATION_STRING "__miss__"

#define MM_OPTIX_DIRECT_CALLABLE_ANNOTATION_STRING "__direct_callable__"

#define MM_OPTIX_CONTINUATION_CALLABLE_ANNOTATION_STRING "__continuation_callable__"

#define MM_OPTIX_EXCEPTION_ANNOTATION_STRING "__exception__"

#define MM_OPTIX_BOUNDS_ANNOTATION_STRING "__boundsKernel__"

namespace megamol {
namespace hpg {
    namespace optix {
        typedef struct box3f {
            glm::vec3 lower;
            glm::vec3 upper;
        } box3f;
    } // namespace optix
} // namespace hpg
} // namespace megamol
