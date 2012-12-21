SET(filter_HEADERS
  filter/ExporterEnsight.hpp
  filter/MeshWriter.hpp
  filter/ExporterHDF5Mesh3D.hpp
  filter/GetPot.hpp
  filter/Exporter.hpp
  filter/ImporterMesh3D.hpp
  filter/Importer.hpp
  filter/ExporterEmpty.hpp
  filter/ExporterVTK.hpp
  filter/ExporterHDF5.hpp
  filter/ImporterMesh2D.hpp
  filter/PartitionIO.hpp
  filter/ParserGmsh.hpp
  filter/ParserINRIAMesh.hpp
CACHE INTERNAL "")

SET(filter_SOURCES
  filter/Importer.cpp
  filter/ImporterMesh3D.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${filter_HEADERS})
  STRING(REPLACE "filter/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/core/filter"
  COMPONENT ${PACKAGE_NAME}
)
