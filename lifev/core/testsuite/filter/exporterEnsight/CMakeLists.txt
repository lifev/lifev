INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ensightExport
  SOURCES ensightExport.cpp ../importExport/RossEthierSteinmanDec.cpp
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  ensightImport
  SOURCES ensightImport.cpp ../importExport/RossEthierSteinmanDec.cpp
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataExportEnsight
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

