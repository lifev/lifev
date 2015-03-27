
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_ventricle
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_ventricle_data
  CREATE_SYMLINK
  SOURCE_FILES MonodomainSolverParamList.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(idealized_test_ventricle
  SOURCE_FILES idealized.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/electrophysiology/data/mesh/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(fiberstest_ventricle
  SOURCE_FILES FiberDirection.h5
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/electrophysiology/data/fibers/
)
