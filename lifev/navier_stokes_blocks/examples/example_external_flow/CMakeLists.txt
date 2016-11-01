
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  external_flow
  SOURCES main.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_external_flow
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(solverOptions_external_flow
  CREATE_SYMLINK
  SOURCE_FILES solversOptionsFast.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_external_flow
  SOURCE_FILES flow_square_coarse.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/navier_stokes_blocks/meshes/
)