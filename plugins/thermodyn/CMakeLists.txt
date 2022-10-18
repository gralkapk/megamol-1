# MegaMol
# Copyright (c) 2021, MegaMol Dev Team
# All rights reserved.
#

megamol_plugin(thermodyn
  BUILD_DEFAULT ON
  DEPENDS_PLUGINS
    mmstd
    geometry_calls
    datatools
    mesh)

if (thermodyn_PLUGIN_ENABLED)
  find_package(PCL CONFIG REQUIRED)
  target_link_libraries(thermodyn PRIVATE ${PCL_LIBRARIES})
endif()
