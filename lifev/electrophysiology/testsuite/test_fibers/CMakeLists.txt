
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_fibers
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_test_fibers_data
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ellipsoidmesh_test_fibers
  SOURCE_FILES idealized.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/electrophysiology/data/mesh/
)

#TRIBITS_COPY_FILES_TO_BINARY_DIR(fibers_test_fibers
#  SOURCE_FILES ellipsoid_fibers
#  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/electrophysiology/data/fibers/
#)
