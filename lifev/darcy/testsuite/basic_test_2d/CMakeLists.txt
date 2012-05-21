
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BasicTest_2d
  SOURCES main.cpp  user_fun.cpp darcy.cpp
  ARGS -c
  NUM_MPI_PROCS 3
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_darcy_basic_test_2d
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)
