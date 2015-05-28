
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  example_benchmarkPurkinje
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_benchmarkPurkinje_data
  CREATE_SYMLINK
  SOURCE_FILES MonodomainSolverParamList.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_benchmarkPurkinje
  SOURCE_FILES benchmark_05mm.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/electrophysiology/data/mesh/
)
