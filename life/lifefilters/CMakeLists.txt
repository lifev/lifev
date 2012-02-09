add_library ( lifefilters ImporterMesh3D.cpp
Exporter.cpp
Importer.cpp )
install ( TARGETS lifefilters DESTINATION lib )
install ( FILES ExporterEmpty.hpp
ExporterEnsight.hpp
ExporterHDF5.hpp
ExporterHDF5Mesh3D.hpp
Exporter.hpp
ExporterVTK.hpp
GetPot.hpp
Importer.hpp
ImporterMesh2D.hpp
ImporterMesh3D.hpp
MeshWriter.hpp DESTINATION include/life/lifefilters)
