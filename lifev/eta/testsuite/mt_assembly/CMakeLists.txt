INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MT_Assembly
  SOURCES main.cpp
  ARGS "10 2"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )
