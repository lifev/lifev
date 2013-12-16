INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  OfflinePartitionIO_Write
  Name OfflinePartitionIO_Write_MeshPartitioner
  SOURCES main_write.cpp
  ARGS "--num-elem 9 --num-parts 3 --partitioner-type MeshPartitioner"
  NUM_MPI_PROCS 1
  COMM mpi
)

TRIBITS_ADD_TEST(
  OfflinePartitionIO_Write
  NAME OfflinePartitionIO_Write_MeshPartitionTool_ParMETIS
  ARGS "--num-elem 9 --num-parts 3 --partitioner-type MeshPartitionTool_ParMETIS"
  NUM_MPI_PROCS 1
  COMM mpi
)

TRIBITS_ADD_TEST(
  OfflinePartitionIO_Write
  NAME OfflinePartitionIO_Write_MeshPartitionTool_Zoltan
  ARGS "--num-elem 9 --num-parts 3 --partitioner-type MeshPartitionTool_Zoltan"
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
