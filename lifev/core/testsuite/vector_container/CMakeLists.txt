INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  VectorContainer
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
  )
