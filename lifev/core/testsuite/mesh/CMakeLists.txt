INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MeshBasic
  SOURCES basic_test.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_meshtest
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_testmesh
  SOURCE_FILES cube_47785.m++
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/mesh++
)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BareMesh
  SOURCES test_baremesh.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_BareMesh
  SOURCE_FILES data_baremesh
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_BareMesh
  SOURCE_FILES cartesian_cube8.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/inria
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  EntitySelection
  SOURCES entity_selection.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_GmshParser
  SOURCE_FILES curve_1d_q1.msh hollow_2d_q2.msh cylinder_3d_p1.msh karman_3d_q2.msh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/gmsh
  CREATE_SYMLINK
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GmshParser
  SOURCES test_gmsh_parser.cpp
  ARGS "curve_1d_q1.msh hollow_2d_q2.msh cylinder_3d_p1.msh karman_3d_q2.msh"
  NUM_MPI_PROCS 1
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_Q2Mesh
  SOURCE_FILES hollow_2d_q1.msh hollow_2d_q2.msh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/gmsh
  CREATE_SYMLINK
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Q2Mesh
  SOURCES test_q2_mesh.cpp
  ARGS "hollow_2d_q1.msh hollow_2d_q2.msh"
  NUM_MPI_PROCS 1
  COMM serial mpi
)

#ADD_SUBDIRECTORY(mesh_partition_tool)
