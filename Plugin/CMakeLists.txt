paraview_add_plugin(LAADFilter
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES LAADFilters
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/LAADFilters/vtk.module")

option(BUILD_TESTING "Build Plugin Testing" OFF)
if (BUILD_TESTING AND BUILD_SHARED_LIBS)
  add_subdirectory(Testing)
endif()

