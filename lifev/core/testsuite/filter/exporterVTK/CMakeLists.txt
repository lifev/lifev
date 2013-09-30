INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  vtkExport
  SOURCES vtkExport.cpp ../importExport/RossEthierSteinmanDec.cpp
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  vtkImport
  SOURCES vtkImport.cpp ../importExport/RossEthierSteinmanDec.cpp
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataExportVTK
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

