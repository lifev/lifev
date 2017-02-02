INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE(
  testExportImport
  SOURCES testExportImport.cpp ../importExport/RossEthierSteinmanDec.cpp
  )

IF (HAVE_HDF5)

TRIBITS_ADD_TEST(
  testExportImport
  POSTFIX_AND_ARGS_0 HDF5 -e hdf5
  NUM_MPI_PROCS 2
  COMM mpi
  )

ENDIF ()


TRIBITS_ADD_TEST(
  testExportImport
  POSTFIX_AND_ARGS_0 Ensight -e ensight -f dataP1
  NUM_MPI_PROCS 2
  COMM mpi
  )

TRIBITS_ADD_TEST(
  testExportImport
  POSTFIX_AND_ARGS_0 VTK -e vtk
  NUM_MPI_PROCS 2
  COMM mpi
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(dataExportEnsight
  CREATE_SYMLINK
  SOURCE_FILES data dataP1
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

