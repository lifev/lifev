INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExporterEnsightToHDF5
  SOURCES main.cpp ensightToHdf5.cpp
  ARGS -c
  NUM_MPI_PROCS 4
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_ExporterEnsightToHDF5
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_ExporterEnsightToHDF5
  SOURCE_FILES cube4x4.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/navier_stokes/data/mesh/inria
)

FILE(GLOB dataFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} importDir/*)

file(MAKE_DIRECTORY importDir)

TRIBITS_COPY_FILES_TO_BINARY_DIR(importDir_ExporterEnsightToHDF5
  SOURCE_FILES "${dataFiles}"
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

