add_library ( lifemesh RegionMesh3DStructured.cpp
ElementShapes.cpp
MeshElementBare.cpp
InternalEntitySelector.cpp
MeshData.cpp
RegionMesh.cpp
MeshVertex.cpp
MeshUtility.cpp
Marker.cpp )
install ( TARGETS lifemesh DESTINATION lib )
install ( FILES BareMesh.hpp
ElementShapes.hpp
InternalEntitySelector.hpp
MarkerDefinitions.hpp
Marker.hpp
MeshChecks.hpp
MeshData.hpp
MeshElementBare.hpp
MeshElement.hpp
MeshElementMarked.hpp
MeshEntity.hpp
MeshPartitioner.hpp
MeshPartitionerOfflineFSI.hpp
MeshUtility.hpp
MeshVertex.hpp
RegionMesh1DBuilders.hpp
RegionMesh3DStructured.hpp
RegionMesh.hpp DESTINATION include/life/lifemesh)
