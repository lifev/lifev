INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tube
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 4
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_tube
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(solverOptions_tube
  CREATE_SYMLINK
  SOURCE_FILES solversOptionsFast.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(solverOptions_tube_rbf
  CREATE_SYMLINK
  SOURCE_FILES SolverParamList_rbf3d.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_fluid_tube
  SOURCE_FILES tube20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_fluid_tube_rotated
  SOURCE_FILES tube20_rotated.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_structure_tube
  SOURCE_FILES vessel20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria
)

