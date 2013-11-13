INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  StaticGraph
  SOURCES main.cpp
  ARGS "10"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )
