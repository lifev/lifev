
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  tube_restart
  SOURCES main.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_tube_restart
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(solverOptions_tube_restart
  CREATE_SYMLINK
  SOURCE_FILES solversOptionsFast.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(solverOptions_tube_rbf_restart
  CREATE_SYMLINK
  SOURCE_FILES SolverParamList_rbf3d.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_fluid_tube_restart
  SOURCE_FILES tube20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_fluid_tube_rotated_restart
  SOURCE_FILES tube20_rotated.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_structure_tube_restart
  SOURCE_FILES vessel20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(fluid_results_h5_restart_P1
CREATE_SYMLINK
SOURCE_FILES fluid_P1P1.h5
SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi_blocks/testsuite/fsi_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(fluid_results_xmf_restart_P1
CREATE_SYMLINK
SOURCE_FILES fluid_P1P1.xmf
SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi_blocks/testsuite/fsi_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(structure_results_h5_restart_P1
CREATE_SYMLINK
SOURCE_FILES structure_P1.h5
SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi_blocks/testsuite/fsi_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(structure_results_xmf_restart_P1
CREATE_SYMLINK
SOURCE_FILES structure_P1.xmf
SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi_blocks/testsuite/fsi_restart
)
 
TRIBITS_COPY_FILES_TO_BINARY_DIR(fluid_results_h5_restart_P2
CREATE_SYMLINK
SOURCE_FILES fluid_P2P1.h5
SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi_blocks/testsuite/fsi_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(fluid_results_xmf_restart_P2
CREATE_SYMLINK
SOURCE_FILES fluid_P2P1.xmf
SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi_blocks/testsuite/fsi_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(structure_results_h5_restart_P2
CREATE_SYMLINK
SOURCE_FILES structure_P2.h5
SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi_blocks/testsuite/fsi_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(structure_results_xmf_restart_P2
CREATE_SYMLINK
SOURCE_FILES structure_P2.xmf
SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi_blocks/testsuite/fsi_restart
)