INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MeshEntityContainer
  SOURCES test_meshentitycontainer.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  VectorSmall
  SOURCES test_vectorsmall.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GhostHandler
  SOURCES test_ghosthandler.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_neighborsRadius
  SOURCES test_neighborsRadius.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_neighborsCircle
  SOURCES test_neighborsCircle.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_ghosthandler
  SOURCE_FILES data_ghosthandler
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_ghosthandler
  SOURCE_FILES mesh_ghosthandler.msh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/freefem
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_neighborsCircle
  SOURCE_FILES data_neighborsCircle
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unstructured_square
  SOURCE_FILES unstructured_square.msh 
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/freefem
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_neighborsRadius
  SOURCE_FILES data_neighborsRadius
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_tube20_long
  SOURCE_FILES tube20_long.mesh 
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/inria
)

