cmake_minimum_required (VERSION 2.6)

project (OpenVolumeMesh)

# add our macro directory to cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

# include our cmake files
include (ACGCommon)

acg_get_version ()

include (ACGOutput)

if (WIN32)
  add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
  )
endif ()

# ========================================================================
# Call the subdirectories with the projects
# ========================================================================

add_subdirectory (src)
add_subdirectory (examples)
add_subdirectory (documentation)

# ========================================================================
# Bundle generation (Targets exist, now configure them)
# ========================================================================

# display results
acg_print_configure_header (OPENVOLUMEMESH "OpenVolumeMesh")
