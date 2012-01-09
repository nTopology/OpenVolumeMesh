include (ACGCommon)

include_directories (
  ..
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# source code directories
set (directories 
  . 
  OpenVolumeMesh
  OpenVolumeMesh/Utils
  OpenHexMesh
  OpenHexMesh/Utils
)

# collect all header and source files
acg_append_files (headers "*.hh" ${directories})
acg_append_files (sources "*.cc" ${directories})

# Disable Library installation when not building OpenVolumeMesh on its own but as part of another project!
if ( NOT ${PROJECT_NAME} MATCHES "OpenVolumeMesh")
  set(ACG_NO_LIBRARY_INSTALL true)
endif()


if (WIN32)
  # OpenVolumeMesh has no dll exports so we have to build a static library on windows
  acg_add_library (OpenVolumeMesh STATIC ${sources} ${headers})
else ()
    acg_add_library (OpenVolumeMesh SHAREDANDSTATIC ${sources} ${headers})
    set_target_properties (OpenVolumeMesh PROPERTIES VERSION ${OPENVOLUMEMESH_VERSION_MAJOR}.${OPENVOLUMEMESH_VERSION_MINOR}
                                               SOVERSION ${OPENVOLUMEMESH_VERSION_MAJOR}.${OPENVOLUMEMESH_VERSION_MINOR} )

endif ()

if ( ${PROJECT_NAME} MATCHES "OpenVolumeMesh")
  # Add core as dependency before fixbundle 
  if ( WIN32 OR APPLE)
    # let bundle generation depend on targets
    add_dependencies (fixbundle OpenVolumeMesh)
  endif ()
endif()


# Install Header Files (Apple)
if ( NOT ACG_PROJECT_MACOS_BUNDLE AND APPLE )
    FILE(GLOB files_install_OVM         "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/*.hh"  "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/*T.cc" )
    FILE(GLOB files_install_OVM_Utils   "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/Utils/*.hh"  "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/Utils/*T.cc" )
    FILE(GLOB files_install_OHM         "${CMAKE_CURRENT_SOURCE_DIR}/OpenHexMesh/*.hh"  "${CMAKE_CURRENT_SOURCE_DIR}/OpenHexMesh/*T.cc" )
    FILE(GLOB files_install_OHM_Utils   "${CMAKE_CURRENT_SOURCE_DIR}/OpenHexMesh/Utils/*.hh"  "${CMAKE_CURRENT_SOURCE_DIR}/OpenHexMesh/Utils/*T.cc" )
    INSTALL(FILES ${files_install_OVM}          DESTINATION include/OpenVolumeMesh/OpenVolumeMesh )
    INSTALL(FILES ${files_install_OVM_Utils}    DESTINATION include/OpenVolumeMesh/OpenVolumeMesh/Utils )
    INSTALL(FILES ${files_install_OHM}          DESTINATION include/OpenVolumeMesh/OpenHexMesh )
    INSTALL(FILES ${files_install_OHM_Utils}    DESTINATION include/OpenVolumeMesh/OpenHexMesh/Utils )
endif()


# Only install if the project name matches OpenVolumeMesh.
if (NOT APPLE AND ${PROJECT_NAME} MATCHES "OpenVolumeMesh")

# Install Header Files)
install(DIRECTORY . 
	DESTINATION include/OpenVolumeMesh
        FILES_MATCHING 
	PATTERN "*.hh"
        PATTERN "CVS" EXCLUDE
        PATTERN ".svn" EXCLUDE
        PATTERN "tmp" EXCLUDE
	PATTERN "Templates" EXCLUDE
        PATTERN "Debian*" EXCLUDE)

#install Template cc files (required by headers)
install(DIRECTORY . 
	DESTINATION include/OpenVolumeMesh
        FILES_MATCHING 
	PATTERN "*T.cc"
        PATTERN "CVS" EXCLUDE
        PATTERN ".svn" EXCLUDE
	PATTERN "tmp" EXCLUDE
	PATTERN "Templates" EXCLUDE
        PATTERN "Debian*" EXCLUDE)

endif ()


