
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_0DFitzHughNagumo
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_0DFitzHughNagumoModel_data
  CREATE_SYMLINK
  SOURCE_FILES FitzHughNagumoParameters.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)


