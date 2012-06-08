INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  vtkExport
  SOURCES vtkExport.cpp ../importExport/RossEthierSteinmanDec.cpp
  ARGS "-p Export"
  NUM_MPI_PROCS 2
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  vtkImport
  SOURCES vtkImport.cpp ../importExport/RossEthierSteinmanDec.cpp
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataExport
  SOURCE_FILES dataExport
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

