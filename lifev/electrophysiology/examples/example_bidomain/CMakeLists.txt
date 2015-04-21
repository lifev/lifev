
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  example_bidomain
  SOURCES main.cpp ElectroETABidomainSolver.hpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_example_bidomain_data
  CREATE_SYMLINK
  SOURCE_FILES BidomainSolverParamList.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(NiedererBM11.mesh_example_bidomain
  SOURCE_FILES benchmark_05mm.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/electrophysiology/data/mesh/
)
