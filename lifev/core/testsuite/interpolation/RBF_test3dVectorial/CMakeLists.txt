INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_rbf3D_vectorial
  SOURCES rbf3d.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_test_rbf3d_vectorial
  SOURCE_FILES data_rbf3d
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverParamList_rbf3d_vectorial
  SOURCE_FILES SolverParamList_rbf3d.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(vessel20_circ_3d_rbf_vectorial.mesh
  SOURCE_FILES vessel20_circ.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(tube20_3d_rbf_vectorial.mesh
  SOURCE_FILES tube20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
)
