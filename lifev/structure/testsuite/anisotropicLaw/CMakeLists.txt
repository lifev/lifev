INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  anisotropicLaw
  NAME SVK
  SOURCES main.cpp ud_functions.cpp
  ARGS "-f data"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataHolzapfel
  SOURCE_FILES dataHolzapfel
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(StructuredCube4_test_structuralsolver.mesh_anisotropicStructure
  SOURCE_FILES StructuredCube4_test_structuralsolver.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/structure/data/mesh/inria/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataParametersHolzapfel
  SOURCE_FILES xmlParameters.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)


