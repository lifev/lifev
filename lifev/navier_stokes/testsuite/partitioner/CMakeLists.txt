INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PartitionerCutAndStore
  SOURCES main_partition.cpp
  ARGS "-f data_partition"
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_partition
  SOURCE_FILES data_partition
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(tube20.mesh_for_partitioner
  SOURCE_FILES tube20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PartitionerLoadAndSolve
  SOURCES main_cylinder.cpp cylinder.cpp
  ARGS "-f data_cylinder"
  NUM_MPI_PROCS 4
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_cylinder
  SOURCE_FILES data_cylinder
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)
