INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Structure
  NAME testStress
  SOURCES main.cpp
  ARGS "-f data"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataStress
  SOURCE_FILES dataSVK
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(StructuredCube4_test_structuralsolver.mesh_Structure
  SOURCE_FILES StructuredCube4_test_structuralsolver.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/structure/data/mesh/inria/
)
