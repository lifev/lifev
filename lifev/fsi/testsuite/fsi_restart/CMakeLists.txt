INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  restartFSI
  SOURCES main.cpp ud_functions.cpp flowConditions.cpp
  NAME restartFSIGI
  ARGS "-f dataCE"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_test_restartFSI
  SOURCE_FILES dataCE
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataGCE_test_restartFSI
  SOURCE_FILES dataGCE
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(vessel20.mesh_test_restartFSI
  SOURCE_FILES vessel20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(tube20.mesh_test_restartFSI
  SOURCE_FILES tube20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
)

TRIBITS_ADD_TEST(
  restartFSI
  NAME restartFSIGCE
  ARGS "-f dataGCE"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileFluidXMFGI
  SOURCE_FILES fluidGI.xmf
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileFluidH5GI
  SOURCE_FILES fluidGI.h5
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileSolidXMFGCE
  SOURCE_FILES solidGCE.xmf
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileSolidH5GCE
  SOURCE_FILES solidGCE.h5
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileFluidXMFGCE
  SOURCE_FILES fluidGCE.xmf
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileFluidH5GCE
  SOURCE_FILES fluidGCE.h5
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileSolidXMFGI
  SOURCE_FILES solidGI.xmf
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileSolidH5GI
  SOURCE_FILES solidGI.h5
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)
