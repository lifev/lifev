INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RegionMarkerID
  SOURCES basic_test.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )
