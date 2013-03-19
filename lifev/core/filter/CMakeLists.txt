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
  filter/ParserGmsh.hpp
  filter/ParserINRIAMesh.hpp
CACHE INTERNAL "")

IF(TPL_ENABLE_HDF5)
    SET(filter_HEADERS
        ${filter_HEADERS}
        filter/HDF5IO.hpp
        filter/PartitionIO.hpp
    CACHE INTERNAL "")
ENDIF()

SET(filter_SOURCES
  filter/Importer.cpp
  filter/ImporterMesh3D.cpp
CACHE INTERNAL "")

IF(TPL_ENABLE_HDF5)
    SET(filter_SOURCES
        ${filter_SOURCES}
        filter/HDF5IO.cpp
    CACHE INTERNAL "")
ENDIF()

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
