INCLUDE(TribitsAddExecutable)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  anisotropicTraction
  SOURCES main.cpp ud_functions.cpp
  COMM serial mpi
  )
