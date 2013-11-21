INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FEFunction
  SOURCES main.cpp fefct.cpp user_fun.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fe_function
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(cartesian_cube8.mesh_fe_function
  SOURCE_FILES cartesian_cube8.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/inria/
)
