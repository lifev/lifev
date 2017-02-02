INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LinearSolverPreconditioner
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 4
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_LinearSolverPreconditioner
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverParamList_xml_LinearSolverPreconditioner
  SOURCE_FILES SolverParamList.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(PrecSolverParamList_xml_LinearSolverPreconditioner
  SOURCE_FILES PrecSolverParamList.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)
