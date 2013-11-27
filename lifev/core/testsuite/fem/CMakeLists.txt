INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CurrentFE
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  QuadRule
  SOURCES test_quadrule.cpp SetOfFun.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Interpolate
  SOURCES test_interpolate.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_Interpolate
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_Interpolate
  SOURCE_FILES hexahexa_4x4x4.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/inria
)

