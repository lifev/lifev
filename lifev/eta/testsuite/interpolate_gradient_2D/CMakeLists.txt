INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ETA_interpolate_gradient_2D
  SOURCES main.cpp ETA_InterpolateGradient2DTest.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  )
