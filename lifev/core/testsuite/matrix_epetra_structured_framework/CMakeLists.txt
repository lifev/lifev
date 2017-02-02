INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixEpetraStructuredFramework
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 4
  COMM serial mpi
  )
