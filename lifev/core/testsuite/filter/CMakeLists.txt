INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

ADD_SUBDIRECTORIES(
  exporterEnsight
  exporterVTK
  exporterHDF5
  getPot
  translator
)

