SET(mesh_HEADERS
  mesh/MeshEntityContainer.hpp
  mesh/MeshData.hpp
  mesh/RegionMesh.hpp
  mesh/Marker.hpp
  mesh/InternalEntitySelector.hpp
  mesh/MeshEntity.hpp
  mesh/MeshVertex.hpp
  mesh/RegionMesh3DStructured.hpp
  mesh/ElementShapes.hpp
  mesh/MeshPartitioner.hpp
  mesh/MarkerDefinitions.hpp
  mesh/MeshUtility.hpp
  mesh/BareMesh.hpp
  mesh/ConvertBareMesh.hpp
  mesh/MeshElement.hpp
  mesh/MeshElementBare.hpp
  mesh/RegionMesh1DBuilders.hpp
  mesh/RegionMesh1DStructured.hpp
  mesh/MeshChecks.hpp
  mesh/MeshElementMarked.hpp
  mesh/GraphCutterBase.hpp
  mesh/GraphCutterZoltan.hpp
  mesh/GraphCutterParMETIS.hpp
  mesh/GraphUtil.hpp
  mesh/MeshPartitionTool.hpp
  mesh/MeshPartBuilder.hpp
  mesh/NeighborMarker.hpp
  mesh/RegionMesh2DStructured.hpp
  mesh/MeshColoring.hpp
CACHE INTERNAL "")

SET(mesh_SOURCES
  mesh/Marker.cpp
  mesh/RegionMesh3DStructured.cpp
  mesh/MeshVertex.cpp
  mesh/RegionMesh.cpp
  mesh/MeshElementBare.cpp
  mesh/ElementShapes.cpp
  mesh/MeshUtility.cpp
  mesh/MeshData.cpp
  mesh/InternalEntitySelector.cpp
  mesh/MeshEntity.cpp
  mesh/RegionMesh2DStructured.cpp
  mesh/MeshColoring.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${mesh_HEADERS})
  STRING(REPLACE "mesh/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/core/mesh"
  COMPONENT ${PACKAGE_NAME}
)
