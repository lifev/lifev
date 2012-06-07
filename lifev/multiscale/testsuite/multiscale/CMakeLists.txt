
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  BasicTest
  NAME BasicTest
  SOURCES main.cpp
  )

TRIBITS_ADD_TEST(
  BasicTest
  POSTFIX_AND_ARGS_0 1D      -s 1 -f MultiscaleDatabase/Test_1D.dat      -c 88.582269095417459
  POSTFIX_AND_ARGS_1 Fluid3D -s 1 -f MultiscaleDatabase/Test_Fluid3D.dat -c 1.8919584275659009
  POSTFIX_AND_ARGS_2 FSI3D   -s 1 -f MultiscaleDatabase/Test_FSI3D.dat   -c 13.147029803470305
  NUM_MPI_PROCS 4
  COMM mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(MS_mesh
  CREATE_SYMLINK
  SOURCE_FILES inria
  DEST_FILES Mesh
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data/mesh
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(MS_functions_database
  CREATE_SYMLINK
  SOURCE_FILES physicalFunctions MultiscaleDatabase
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data
)

