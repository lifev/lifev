SET(filter_HEADERS
filter/ExporterEmpty.hpp
filter/ExporterEnsight.hpp
filter/ExporterHDF5.hpp
filter/ExporterHDF5Mesh3D.hpp
filter/Exporter.hpp
filter/ExporterVTK.hpp
filter/GetPot.hpp
filter/Importer.hpp
filter/ImporterMesh2D.hpp
filter/ImporterMesh3D.hpp
filter/MeshWriter.hpp
CACHE INTERNAL "")

SET(filter_SOURCES
filter/Exporter.cpp
filter/Importer.cpp
filter/ImporterMesh3D.cpp
CACHE INTERNAL "")
