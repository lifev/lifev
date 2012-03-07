INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
BareMesh.hpp
ElementShapes.hpp
InternalEntitySelector.hpp
MarkerDefinitions.hpp
Marker.hpp
MeshChecks.hpp
MeshData.hpp
MeshElementBare.hpp
MeshElement.hpp
MeshElementMarked.hpp
MeshEntityContainer.hpp
MeshEntity.hpp
MeshPartitioner.hpp
MeshUtility.hpp
MeshVertex.hpp
RegionMesh1DBuilders.hpp
RegionMesh3DStructured.hpp
RegionMesh.hpp
)

SET(SOURCES ${SOURCES}
ElementShapes.cpp
InternalEntitySelector.cpp
Marker.cpp
MeshData.cpp
MeshElementBare.cpp
MeshUtility.cpp
MeshVertex.cpp
RegionMesh3DStructured.cpp
RegionMesh.cpp
)
