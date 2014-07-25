
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_benchmark
  SOURCES main.cpp benchmarkUtility.hpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(data_benchmark_data
  CREATE_SYMLINK
  SOURCE_FILES MonodomainSolverParamList.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_benchmark
  SOURCE_FILES benchmark_05mm.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/electrophysiology/data/mesh/
)
