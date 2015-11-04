INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

IF(LevelSet_ENABLE_ETA)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  1_Volume_computation
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

ENDIF()
