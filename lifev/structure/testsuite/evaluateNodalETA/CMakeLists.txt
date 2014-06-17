INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  evaluateNodalETA
  NAME evaluateNodal
  SOURCES main.cpp ud_functions.cpp
  ARGS "-f data"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataEvaluate
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(cubeY4.mesh_evaluateNodal
  SOURCE_FILES cubeY4.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/structure/data/mesh/inria/
)
