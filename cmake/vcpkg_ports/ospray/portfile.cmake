set(OSPRAY_VERSION 2.10.0)

vcpkg_from_github(
  OUT_SOURCE_PATH SOURCE_PATH_OSPRAY
  REPO ospray/OSPRay
  REF v${OSPRAY_VERSION}
  SHA512 b72e6dfc6bd8c7f4013e7074d89a47d52ff56e1270828b8cafc7365209c68614ce4e5842218ff30b170e2ab206741f9c57bf3d31649c5e70785a01267d073e62
  HEAD_REF master
  PATCHES
    fix-config-dir-depth.patch
)

vcpkg_from_github(
  OUT_SOURCE_PATH SOURCE_PATH_PKD
  REPO UniStuttgart-VISUS/ospray-module-pkd
  REF 8694e0de85a3d1719f1dc062f7f41967d6d8e408
  SHA512 1e14938aed9ba939cb0b0e9c4af0787713774b914448ced9e66d476d86b8fec0b052b3ad2c2f74fc61ec5f272ec4d8a0ee44da70b573e7dbccca2f7608e7baee
  HEAD_REF megamol_osp2
)

file(COPY "${SOURCE_PATH_PKD}" DESTINATION "${SOURCE_PATH_OSPRAY}/modules")

# currently this does not work, either you decide for one or you take them all
# if (sse4 IN_LIST FEATURES)
#   list(APPEND OSPRAY_ISA SSE4)
# endif()
# if(avx IN_LIST FEATURES)
#   list(APPEND OSPRAY_ISA AVX)
# endif()
# if(avx2 IN_LIST FEATURES)
#   list(APPEND OSPRAY_ISA AVX2)
# endif()
# if(avx512 IN_LIST FEATURES)
#   list(APPEND OSPRAY_ISA AVX512)
# endif()

vcpkg_cmake_configure(
  SOURCE_PATH "${SOURCE_PATH_OSPRAY}"
  OPTIONS
    -DOSPRAY_BUILD_ISA=ALL
    -DOSPRAY_ENABLE_APPS=false
    -DOSPRAY_MODULE_PKD=ON
)
vcpkg_cmake_install()

vcpkg_copy_pdbs()
vcpkg_cmake_config_fixup()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/share")

# file(
#   INSTALL "${SOURCE_PATH}/LICENSE"
#   DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}"
#   RENAME copyright)

file(INSTALL "${SOURCE_PATH_OSPRAY}/LICENSE.txt" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}" RENAME copyright)
