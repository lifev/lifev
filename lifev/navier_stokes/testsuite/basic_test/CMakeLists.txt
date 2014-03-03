
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BasicTest
  SOURCES mainRossEthierSteinman.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataRossEthierSteinman
  SOURCE_FILES dataRossEthierSteinman
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(cube4x4.mesh
  SOURCE_FILES cube4x4.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/navier_stokes/data/mesh/inria/
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BasicTest2D
  SOURCES mainKimMoin.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataKimMoin
  SOURCE_FILES dataKimMoin
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(square20x20.msh
  SOURCE_FILES square20x20.msh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/freefem/
)

