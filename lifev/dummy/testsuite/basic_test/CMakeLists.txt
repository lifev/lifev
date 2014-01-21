
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  basic_test
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_COPY_FILES_TO_BINARY_DIR(dataDummy
#   SOURCE_FILES data
#   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
# )

# TRIBITS_COPY_FILES_TO_BINARY_DIR(cartesian_cube8.mesh
#   SOURCE_FILES cartesian_cube8.mesh
#   SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/darcy/data/mesh/inria
# )

