INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TemplateTest
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_template_test
  SOURCE_FILES data.txt
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

# The previous command can be used to copy a mesh for a test.
# For example:
# TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_template_test
#   SOURCE_FILES some_mesh_name.mesh
#   SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/inria
# )
