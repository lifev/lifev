
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  benchmark_GW
  SOURCES main.cpp
  COMM serial mpi
  )
