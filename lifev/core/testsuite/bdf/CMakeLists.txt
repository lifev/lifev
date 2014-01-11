INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BDF2D
  SOURCES main.cpp test_bdf_2d.cpp
  NAME BDF2D-Order1
  ARGS "-f data_2d_order1"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_bdf_2d_order1
  SOURCE_FILES data_2d_order1
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(bdf_2d_mesh
  SOURCE_FILES square20x20.msh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/freefem/
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BDF3D
  SOURCES main.cpp test_bdf.cpp
  NAME BDF3D-Order2
  ARGS "-f data_3d_order2"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_bdf_3d_order2
  SOURCE_FILES data_3d_order2
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_ADD_TEST(
  BDF2D
  NAME BDF2D-Order3
  ARGS "-f data_2d_order3"
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_bdf_2d_order3
  SOURCE_FILES data_2d_order3
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(bdf_3d_mesh
  SOURCE_FILES cube4x4.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/inria
)
