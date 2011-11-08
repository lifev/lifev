add_library (lifemesh ElementShapes.cpp InternalEntitySelector.cpp Marker.cpp MeshData.cpp MeshElementBare.cpp MeshUtility.cpp MeshVertex.cpp RegionMesh3D.cpp RegionMesh3DStructured.cpp)
install (TARGETS lifemesh DESTINATION lib)
install (FILES BareMesh.hpp ElementShapes.hpp InternalEntitySelector.hpp MarkerDefinitions.hpp Marker.hpp MeshChecks.hpp MeshData.hpp MeshElementBare.hpp MeshElement.hpp MeshElementMarked.hpp MeshEntity.hpp MeshPartitioner.hpp MeshPartitionerOfflineFSI.hpp MeshUtility.hpp MeshVertex.hpp RegionMesh1D.hpp RegionMesh2D.hpp RegionMesh3D.hpp RegionMesh3DStructured.hpp DESTINATION include/life/lifemesh)

