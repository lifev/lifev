INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MeshPartitionTool
  NAME MeshPartitionTool_ParMETIS
  SOURCES main.cpp
  ARGS "--num-elem 9 --graph-lib parmetis"
  NUM_MPI_PROCS 3
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MeshPartitionTool
  NAME MeshPartitionTool_Zoltan
  ARGS "--num-elem 9 --graph-lib zoltan"
  NUM_MPI_PROCS 3
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
