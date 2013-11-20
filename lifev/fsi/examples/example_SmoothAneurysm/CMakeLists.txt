
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  SmoothAneurysm
  SOURCES main.cpp ud_functions.cpp lumpedHeart.cpp flowConditions.cpp
  COMM serial mpi
  )

# TRIBITS_COPY_FILES_TO_BINARY_DIR(data_exampleSmoothAneurysm
#   SOURCE_FILES data
#   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
#   CREATE_SYMLINK
# )

# TRIBITS_COPY_FILES_TO_BINARY_DIR(fluidBig.mesh_test_monolitic
#   SOURCE_FILES fluidBig.mesh
# #  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
#   CREATE_SYMLINK
# )

# TRIBITS_COPY_FILES_TO_BINARY_DIR(solidBig.mesh_test_monolitic
#   SOURCE_FILES solidBig.mesh
# #  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
#   CREATE_SYMLINK
# )
