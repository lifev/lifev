add_library ( lifefilters Exporter.cpp
Importer.cpp
ImporterMesh3D.cpp )
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
