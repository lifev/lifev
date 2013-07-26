INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  anisotropicLaw
  NAME Holzapfel
  SOURCES main.cpp ud_functions.cpp
  ARGS "-f dataHolzapfel"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataHolzapfel
  SOURCE_FILES dataHolzapfel
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(cubeY4.mesh_anisotropicStructure
  SOURCE_FILES cubeY4.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/structure/data/mesh/inria/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataParametersHolzapfel
  SOURCE_FILES xmlParameters.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_ADD_TEST(
  anisotropicLaw
  NAME Multimechanism
  ARGS "-f dataMultimechanism"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataMultimechanismTest
  SOURCE_FILES dataMultimechanism
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)

TRIBITS_ADD_TEST(
  anisotropicLaw
  NAME DistributedAnisotropic
  ARGS "-f dataDistributed"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataDistributedTest
  SOURCE_FILES dataDistributed
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)


