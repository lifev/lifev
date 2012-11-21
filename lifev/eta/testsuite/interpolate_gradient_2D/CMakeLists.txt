INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ETA_VectorialADR_2D
  SOURCES main.cpp ETA_VectorialADR2DTest.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  )
