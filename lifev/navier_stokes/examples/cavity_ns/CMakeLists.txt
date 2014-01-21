
INCLUDE(TribitsAddExecutable)

# Add executable called "cavity_ns_test" that is built from the source files
# "cavity_ns.cpp".
TRIBITS_ADD_EXECUTABLE(
  cavity_ns_test
  SOURCES cavity_ns.cpp
  COMM serial mpi
  )

# Access the data
TRIBITS_COPY_FILES_TO_BINARY_DIR(
  data_cavity_ns_test
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

# Access the mesh
TRIBITS_COPY_FILES_TO_BINARY_DIR(
  meshes_cavity_ns_test
  SOURCE_FILES cube4x4.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/freefem/
)
