INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  3_ETA_rhs_and_value
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  )
