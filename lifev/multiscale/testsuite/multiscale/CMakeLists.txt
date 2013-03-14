
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  Framework
  NAME Framework
  SOURCES main.cpp
  )

IF (LifeV_ENABLE_NavierStokes)

TRIBITS_ADD_TEST(
  Framework
  POSTFIX_AND_ARGS_0 Fluid3D_MTNS       -s 1 -f MultiscaleDatabase/Test_Fluid3D.dat                   -c 267617.11076637130 -o Test_1
  NUM_MPI_PROCS 4
  COMM mpi
  )
  
ENDIF ()

IF (LifeV_ENABLE_FSI)

TRIBITS_ADD_TEST(
  Framework
  POSTFIX_AND_ARGS_0 FSI3D_MTNS         -s 1 -f MultiscaleDatabase/Test_FSI3D.dat                     -c 69.503547600977484 -o Test_3
  NUM_MPI_PROCS 4
  COMM mpi
  )
  
ENDIF ()

IF (LifeV_ENABLE_OneDFSI)
  
TRIBITS_ADD_TEST(
  Framework
  POSTFIX_AND_ARGS_0 FSI1D_MNS          -s 1 -f MultiscaleDatabase/Test_FSI1D.dat                     -c 1281.4882585623886 -o Test_2
  NUM_MPI_PROCS 4
  COMM mpi
  )
  
IF (LifeV_ENABLE_ZeroDimensional)
TRIBITS_ADD_TEST(
  Framework
  POSTFIX_AND_ARGS_0 Windkessel0D       -s 1 -f MultiscaleDatabase/Test_Windkessel0D.dat              -c 105330.10083689842 -o Test_0
  NUM_MPI_PROCS 4
  COMM mpi
  )
ENDIF ()

IF (LifeV_ENABLE_FSI)
TRIBITS_ADD_TEST(
  Framework
  POSTFIX_AND_ARGS_0 1D-FSI3D-1D_MNSA   -s 1 -f MultiscaleDatabase/Test_FSI1D-FSI3D-FSI1D.dat         -c 1278.3958233134831 -o Test_4
  POSTFIX_AND_ARGS_1 Restart            -s 1 -f MultiscaleDatabase/Test_FSI1D-FSI3D-FSI1D_Restart.dat -c 506.64046128381432 -o Test_4
  NUM_MPI_PROCS 4
  COMM mpi
  )
ENDIF ()
  
ENDIF ()


TRIBITS_COPY_FILES_TO_BINARY_DIR(MS_mesh
  CREATE_SYMLINK
  SOURCE_FILES inria
  DEST_FILES Mesh
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data/mesh
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(MS_functions_database
  CREATE_SYMLINK
  SOURCE_FILES FunctionsDatabase MultiscaleDatabase
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data
)