# MegaMol
# Copyright (c) 2021, MegaMol Dev Team
# All rights reserved.
#

megamol_plugin(optix_hpg
  BUILD_DEFAULT OFF
  DEPENDS_FEATURES
    cuda
  DEPENDS_PLUGINS
    mmstd
    mesh
    geometry_calls)

if (optix_hpg_PLUGIN_ENABLED)
  find_package(nvcomp CONFIG REQUIRED)
  find_package(glm CONFIG REQUIRED)
  target_link_libraries(optix_hpg
    PRIVATE
      glm::glm
      nvcomp::nvcomp)

  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

  set(OPTIX_HPG_CUDA_ARCH "75" CACHE STRING "CUDA Arch for Optix HPG")
  set_property(CACHE OPTIX_HPG_CUDA_ARCH PROPERTY STRINGS 61 75 80)

  if (OPTIX_HPG_CUDA_ARCH STREQUAL "61")
    set(CMAKE_CUDA_ARCHITECTURES 61)
  elseif (OPTIX_HPG_CUDA_ARCH STREQUAL "75")
    set(CMAKE_CUDA_ARCHITECTURES 75)
  elseif (OPTIX_HPG_CUDA_ARCH STREQUAL "80")
    set(CMAKE_CUDA_ARCHITECTURES 80)
  endif ()

  set(OPTIX_POWER_DEC_VAL "4" CACHE STRING "Fixed Point Prec")
  target_compile_definitions(optix_hpg PRIVATE MM_OPTIX_POWER_DEC_VAL=${OPTIX_POWER_DEC_VAL})
  add_definitions(-DMM_OPTIX_POWER_DEC_VAL=${OPTIX_POWER_DEC_VAL})

  include(configure_optix)
  find_package(CUDAToolkit)

  get_target_property(GLM_INCLUDES glm::glm INTERFACE_INCLUDE_DIRECTORIES)

  # Collect source files
  file(GLOB_RECURSE cuda_resource_include_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "cuda_resources/*.h" "cuda_resources/*.cuh")

  include_directories(${GLM_INCLUDES} "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/cuda_resources")
  cuda_compile_and_embed(embedded_raygen_programs cuda_resources/raygen.cu)
  cuda_compile_and_embed(embedded_miss_programs cuda_resources/miss.cu)
  cuda_compile_and_embed(embedded_miss_occlusion_programs cuda_resources/miss_occlusion.cu)
  cuda_compile_and_embed(embedded_sphere_programs cuda_resources/sphere.cu)
  cuda_compile_and_embed(embedded_sphere_occlusion_programs cuda_resources/sphere_occlusion.cu)
  cuda_compile_and_embed(embedded_mesh_programs cuda_resources/mesh.cu)
  cuda_compile_and_embed(embedded_transitioncalculator_programs cuda_resources/transitioncalculator.cu)
  cuda_compile_and_embed(embedded_pkd_programs cuda_resources/pkd.cu)

  target_sources(optix_hpg PRIVATE ${cuda_resource_include_files} ${embedded_raygen_programs} ${embedded_miss_programs} ${embedded_miss_occlusion_programs} ${embedded_sphere_programs} ${embedded_sphere_occlusion_programs} ${embedded_mesh_programs} ${embedded_transitioncalculator_programs} ${embedded_pkd_programs})

  target_include_directories(optix_hpg PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/cuda_resources/")
  target_link_libraries(optix_hpg PRIVATE CUDA::cudart CUDA::cuda_driver)

  if (WIN32)
    get_target_property(NVCOMP_DLL nvcomp::nvcomp IMPORTED_LOCATION_RELEASE)
    message(INFO " ${NVCOMP_DLL}")
    if (NVCOMP_DLL MATCHES ".dll$")
      get_filename_component(NVCOMP_DLL_PATH ${NVCOMP_DLL} DIRECTORY)
      message(INFO " ${NVCOMP_DLL_PATH}")
      file(GLOB NVCOMP_DLLS "${NVCOMP_DLL_PATH}/*.dll")
      message(INFO " ${NVCOMP_DLLS}")
      install(FILES ${NVCOMP_DLLS} DESTINATION bin)
    endif ()
  endif ()
  if (MEGAMOL_USE_OPENGL)
    target_link_libraries(optix_hpg PRIVATE mmstd_gl)
  endif ()
endif ()
