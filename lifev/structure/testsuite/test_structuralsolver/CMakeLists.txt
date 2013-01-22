INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Structure
  NAME SVK
  SOURCES main.cpp
  ARGS "-f dataSVK"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataStVenant
  SOURCE_FILES dataSVK
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(StructuredCube4_test_structuralsolver.mesh_Structure
  SOURCE_FILES StructuredCube4_test_structuralsolver.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/structure/data/mesh/inria/
)

TRIBITS_ADD_TEST(
  Structure
  NAME NH
  ARGS "-f dataNH"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataNeoHookean
  SOURCE_FILES dataNH
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_ADD_TEST(
  Structure
  NAME LE
  ARGS "-f dataLE"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataLinearElastic
  SOURCE_FILES dataLE
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_ADD_TEST(
  Structure
  NAME EXP
  ARGS "-f dataEXP"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataExponential
  SOURCE_FILES dataEXP
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)
