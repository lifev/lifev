SET(mesh_HEADERS
mesh/BareMesh.hpp
mesh/ElementShapes.hpp
mesh/InternalEntitySelector.hpp
mesh/MarkerDefinitions.hpp
mesh/Marker.hpp
mesh/MeshChecks.hpp
mesh/MeshData.hpp
mesh/MeshElementBare.hpp
mesh/MeshElement.hpp
mesh/MeshElementMarked.hpp
mesh/MeshEntityContainer.hpp
mesh/MeshEntity.hpp
mesh/MeshPartitioner.hpp
mesh/MeshUtility.hpp
mesh/MeshVertex.hpp
mesh/RegionMesh1DBuilders.hpp
mesh/RegionMesh3DStructured.hpp
mesh/RegionMesh.hpp
CACHE INTERNAL "")

SET(mesh_SOURCES
mesh/ElementShapes.cpp
mesh/InternalEntitySelector.cpp
mesh/Marker.cpp
mesh/MeshData.cpp
mesh/MeshElementBare.cpp
mesh/MeshUtility.cpp
mesh/MeshVertex.cpp
mesh/RegionMesh3DStructured.cpp
mesh/RegionMesh.cpp
CACHE INTERNAL "")
