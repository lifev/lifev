INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TensionsDispl
  NAME StressDispl
  SOURCES main.cpp
  ARGS "-f dataDisplacement"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataStress
  SOURCE_FILES dataDisplacement
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(cubeY4.mesh_Structure
  SOURCE_FILES cubeY4.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/structure/data/mesh/inria/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(xmfDataStructure
  SOURCE_FILES structure.xmf
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/structure/testsuite/test_evaluationStresses/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(h5DataStructure
  SOURCE_FILES structure.h5
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/structure/testsuite/test_evaluationStresses/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataTensions
  SOURCE_FILES dataTensions
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataEigenvalues
  SOURCE_FILES dataEigenvalues
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)
