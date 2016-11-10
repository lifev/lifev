INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  nnlin_time_dep_ns
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_nnlin_time_dep_ns
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_nnlin_time_dep_ns
SOURCE_FILES tube20.mesh
SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(solverOptions_nnlin_time_dep_ns
  CREATE_SYMLINK
  SOURCE_FILES solversOptionsFast.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

