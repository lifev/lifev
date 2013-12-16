INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PartitionIO_Write
  SOURCES main_write.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_partition_io
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PartitionIO_Read
  SOURCES main_read.cpp
  ARGS -c
  NUM_MPI_PROCS 3
  COMM mpi
)
