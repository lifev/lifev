INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  restartFSI
  SOURCES main.cpp ud_functions.cpp flowConditions.cpp
  NAME restartFSI
  ARGS "-f dataCE"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_test_restartFSI
  SOURCE_FILES dataCE
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(vessel20.mesh_test_restartFSI
  SOURCE_FILES vessel20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(tube20.mesh_test_restartFSI
  SOURCE_FILES tube20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileFluidXMF
  SOURCE_FILES fluid.xmf
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileFluidH5
  SOURCE_FILES fluid.h5
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileSolidXMF
  SOURCE_FILES solid.xmf
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataFileSolidH5
  SOURCE_FILES solid.h5
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/testsuite/fsi_restart/
)


