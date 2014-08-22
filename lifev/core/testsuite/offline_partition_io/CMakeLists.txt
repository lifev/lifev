INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  OfflinePartitionIO_Write
  Name OfflinePartitionIO_Write_MeshPartitioner
  SOURCES main_write.cpp
  NUM_MPI_PROCS 1
  COMM mpi
)

TRIBITS_ADD_TEST(
  OfflinePartitionIO_Write
  NAME OfflinePartitionIO_Write_MeshPartitionTool_ParMETIS
  NUM_MPI_PROCS 1
  COMM mpi
)

TRIBITS_ADD_TEST(
  OfflinePartitionIO_Write
  NAME OfflinePartitionIO_Write_MeshPartitionTool_Zoltan
  NUM_MPI_PROCS 1
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  OfflinePartitionIO_Read
  NAME OfflinePartitionIO_Read_MeshPartitioner
  SOURCES main_read.cpp
  ARGS "--partitioner-type MeshPartitioner"
  NUM_MPI_PROCS 3
  COMM mpi
)

TRIBITS_ADD_TEST(
  OfflinePartitionIO_Read
  NAME OfflinePartitionIO_Read_MeshPartitionTool_ParMETIS
  ARGS "--partitioner-type MeshPartitionTool_ParMETIS"
  NUM_MPI_PROCS 3
  COMM mpi
)

TRIBITS_ADD_TEST(
  OfflinePartitionIO_Read
  NAME OfflinePartitionIO_Read_MeshPartitionTool_Zoltan
  ARGS "--partitioner-type MeshPartitionTool_Zoltan"
  NUM_MPI_PROCS 3
  COMM mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_offlineMeshPartitioner
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_offlinePartitioner
  SOURCE_FILES tube20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
)