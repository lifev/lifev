
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  test_block_solver
  SOURCES main.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_test_block_solver
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_test_block_solver
  SOURCE_FILES tube20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria
)

