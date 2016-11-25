#openvolumemesh

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4244")

list(APPEND sources 
    "../src/OpenVolumeMesh/Core/BaseEntities.cc"
    "../src/OpenVolumeMesh/Core/BaseEntities.hh"
    "../src/OpenVolumeMesh/Core/BaseProperty.cc"
    "../src/OpenVolumeMesh/Core/BaseProperty.hh"
    "../src/OpenVolumeMesh/Core/Iterators.cc"
    "../src/OpenVolumeMesh/Core/Iterators.hh"
    "../src/OpenVolumeMesh/Core/OpenVolumeMeshBaseProperty.cc"
    "../src/OpenVolumeMesh/Core/OpenVolumeMeshBaseProperty.hh"
    "../src/OpenVolumeMesh/Core/OpenVolumeMeshHandle.cc"
    "../src/OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh"
    "../src/OpenVolumeMesh/Core/PropertyDefines.cc"
    "../src/OpenVolumeMesh/Core/PropertyDefines.hh"
    "../src/OpenVolumeMesh/Core/ResourceManager.cc"
    "../src/OpenVolumeMesh/Core/ResourceManager.hh"
    "../src/OpenVolumeMesh/Core/Serializers.cc"
    "../src/OpenVolumeMesh/Core/Serializers.hh"
    "../src/OpenVolumeMesh/Core/TopologyKernel.cc"
    "../src/OpenVolumeMesh/Core/TopologyKernel.hh"
    "../src/OpenVolumeMesh/FileManager/FileManager.cc"
    "../src/OpenVolumeMesh/FileManager/FileManager.hh"
    "../src/OpenVolumeMesh/Geometry/Vector.cc"
    "../src/OpenVolumeMesh/Mesh/HexahedralMeshIterators.cc"
    "../src/OpenVolumeMesh/Mesh/HexahedralMeshIterators.hh"
    "../src/OpenVolumeMesh/Mesh/HexahedralMeshTopologyKernel.cc"
    "../src/OpenVolumeMesh/Mesh/HexahedralMeshTopologyKernel.hh"
    "../src/OpenVolumeMesh/Mesh/TetrahedralMeshIterators.cc"
    "../src/OpenVolumeMesh/Mesh/TetrahedralMeshIterators.hh"
    "../src/OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.cc"
    "../src/OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh"
    "../src/OpenVolumeMesh/Attribs/OpenVolumeMeshStatus.cc"
    "../src/OpenVolumeMesh/Attribs/OpenVolumeMeshStatus.hh"
    "../src/OpenVolumeMesh/Attribs/StatusAttrib.cc"
    "../src/OpenVolumeMesh/Attribs/StatusAttrib.hh"
)
add_library(openvolumemesh STATIC ${sources})

#target_link_libraries(openvolumemesh)

target_include_directories(openvolumemesh PUBLIC
    ../src
)

target_compile_definitions(openvolumemesh PUBLIC
    _USE_MATH_DEFINES
    INCLUDE_TEMPLATES
)
