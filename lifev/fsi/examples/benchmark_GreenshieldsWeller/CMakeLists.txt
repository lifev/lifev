
INCLUDE(TribitsAddExecutable)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  benchmark_GW
  SOURCES main.cpp
  COMM serial mpi
  )
