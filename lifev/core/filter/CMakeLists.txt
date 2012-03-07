INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
ExporterEmpty.hpp
ExporterEnsight.hpp
ExporterHDF5.hpp
ExporterHDF5Mesh3D.hpp
Exporter.hpp
ExporterVTK.hpp
GetPot.hpp
Importer.hpp
ImporterMesh2D.hpp
ImporterMesh3D.hpp
MeshWriter.hpp
)

SET(SOURCES ${SOURCES}
Exporter.cpp
Importer.cpp
ImporterMesh3D.cpp
)
