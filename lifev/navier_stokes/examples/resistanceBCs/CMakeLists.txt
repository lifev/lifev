
INCLUDE(TribitsAddExecutable)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  ResistanceBCs
  SOURCES main.cpp cylinder.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataResistance
  CREATE_SYMLINK
  SOURCE_FILES dataResistance
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

# TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_Cylinder
#   SOURCE_FILES tube20.mesh
#   SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
# )

