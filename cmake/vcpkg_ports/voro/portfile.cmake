vcpkg_check_linkage(ONLY_STATIC_LIBRARY)

vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO chr1shr/voro
    REF 56d619faf3479313399516ad71c32773c29be859
    SHA512 e2f223b7c79c3a6dd2d253727d5a60b92d7fde3b0637dfa25657050f558b9dcd6f35d4fd060f70995a4dfeae8a918f7708c3c2daeec6462390f23dd9b4225dd1
    HEAD_REF master
)

# Hotfix the original voro++ CMakeLists.txt to include ALL headers
file(READ "${SOURCE_PATH}/CMakeLists.txt" filedata)
string(REPLACE "src/voro++.hh" "src/*.hh" filedata "${filedata}")
file(WRITE "${SOURCE_PATH}/CMakeLists.txt" "${filedata}")

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DVORO_BUILD_EXAMPLES:BOOL=OFF
        -DVORO_BUILD_CMD_LINE:BOOL=OFF
        -DVORO_ENABLE_DOXYGEN:BOOL=OFF
        -DCMAKE_INSTALL_INCLUDEDIR="${CURRENT_PACKAGES_DIR}/include/${PORT}"
)
vcpkg_cmake_install()
vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/VORO)

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/share")

# Copy the copyright notice
file(INSTALL "${SOURCE_PATH}/LICENSE" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}" RENAME copyright)