INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  11_TestingLaplacianPhiI
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  )
