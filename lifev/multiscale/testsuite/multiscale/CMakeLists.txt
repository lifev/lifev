
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  BasicTest
  NAME BasicTest
  SOURCES main.cpp
  )

TRIBITS_ADD_TEST(
  BasicTest
  POSTFIX_AND_ARGS_0 Fluid3D -s 1 -f MultiscaleDatabase/Test_Fluid3D.dat -c 1.8919584275659009 -o test_Fluid3D
  POSTFIX_AND_ARGS_1 1D      -s 1 -f MultiscaleDatabase/Test_1D.dat      -c 88.582269095417459 -o test_1D
  POSTFIX_AND_ARGS_2 FSI3D   -s 1 -f MultiscaleDatabase/Test_FSI3D.dat   -c 4.3399328834457638 -o test_FSI3D
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

