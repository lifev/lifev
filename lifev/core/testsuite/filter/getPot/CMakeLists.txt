INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GetPot
  SOURCES main_getpot.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_getpot
  SOURCE_FILES data_getpot
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

