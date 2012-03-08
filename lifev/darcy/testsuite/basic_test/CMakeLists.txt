
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  basic_test
  SOURCES main.cpp  user_fun.cpp darcy.cpp
  ARGS -c
  NUM_MPI_PROCS 3
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(cartesian_cube8.mesh
  SOURCE_FILES cartesian_cube8.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/darcy/data/mesh/inria
)

