
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Monolithic
  SOURCES main.cpp ud_functions.cpp lumpedHeart.cpp flowConditions.cpp
  NAME MonolithicGI
  ARGS "-f dataCE"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_test_monolitic_gi
  SOURCE_FILES dataCE
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(vessel20.mesh_test_monolitic
  SOURCE_FILES vessel20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(tube20.mesh_test_monolitic
  SOURCE_FILES tube20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
)

TRIBITS_ADD_TEST(
  Monolithic
  NAME MonolithicGE
  ARGS "-f dataGCE"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_test_monolitic_ge
  SOURCE_FILES dataGCE
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

