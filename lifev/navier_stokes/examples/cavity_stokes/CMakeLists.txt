
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

# Make sure the compiler can find include files
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# Add executable called "cavity_stokes_test" that is built from the source files
# "cavity_stokes.cpp". (Note that the name of the executable has to be uniquely chosen.) 
TRIBITS_ADD_EXECUTABLE(
  cavity_stokes_test
  SOURCES cavity_stokes.cpp
  COMM serial mpi
  )

# Access the data
TRIBITS_COPY_FILES_TO_BINARY_DIR(
  data_cavity_stokes_test
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

# Access the meshes
TRIBITS_COPY_FILES_TO_BINARY_DIR(
  meshes_cavity_stokes_test
  SOURCE_FILES cube4x4.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/freefem/
)
