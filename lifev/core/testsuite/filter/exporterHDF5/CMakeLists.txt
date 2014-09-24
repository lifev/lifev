INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  hdf5Export
  SOURCES hdf5Export.cpp ../importExport/RossEthierSteinmanDec.cpp
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  hdf5Import
  SOURCES hdf5Import.cpp ../importExport/RossEthierSteinmanDec.cpp
  NUM_MPI_PROCS 2
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataExportHdf5
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

