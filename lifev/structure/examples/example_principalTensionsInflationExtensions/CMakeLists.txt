INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TensionsDisplCyl
  NAME StressDisplCyl
  SOURCES main.cpp
  ARGS "-f dataCylindrical"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataCylindrical
  SOURCE_FILES dataCylindrical
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)

