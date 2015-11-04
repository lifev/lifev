INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tensions
  NAME reconstructDisplacement
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
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/structure/testsuite/principalTensions/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(h5DataStructure
  SOURCE_FILES structure.h5
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/structure/testsuite/principalTensions/
)

TRIBITS_ADD_TEST(
  Tensions
  NAME reconstructTensions
  ARGS "-f dataTensions"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataTensions
  SOURCE_FILES dataTensions
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)

TRIBITS_ADD_TEST(
  Tensions
  NAME reconstructEigenvalues
  ARGS "-f dataEigenvalues"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataEigenvalues
  SOURCE_FILES dataEigenvalues
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)

TRIBITS_ADD_TEST(
  Tensions
  NAME VonMises-LE
  ARGS "-f dataLE"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataStressLE
  SOURCE_FILES dataLE
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_ADD_TEST(
  Tensions
  NAME VonMises-SVK
  ARGS "-f dataSVK"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataStressSVK
  SOURCE_FILES dataSVK
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_ADD_TEST(
  Tensions
  NAME VonMises-SVKPenalized
  ARGS "-f dataSVKPenalized"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataStressSVKPenalized
  SOURCE_FILES dataSVKPenalized
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_ADD_TEST(
  Tensions
  NAME VonMises-NH
  ARGS "-f dataNH"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataStressNH
  SOURCE_FILES dataNH
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_ADD_TEST(
  Tensions
  NAME VonMises-2ndExponential
  ARGS "-f data2ndExpo"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataStress2ndExpo
  SOURCE_FILES data2ndExpo
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TensionsAnistropic
  NAME Anisotropic
  SOURCES mainAnisotropic.cpp ud_functions.cpp
  ARGS "-f dataHolzapfel"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataHolzapfelEvaluateStress
  SOURCE_FILES dataHolzapfel
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)
