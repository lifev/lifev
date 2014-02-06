
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  externalFluidFlow
  SOURCES main.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_externalFluidFlow
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_externalFluidFlow
  SOURCE_FILES cylinder.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/eta/testsuite/VMS_Cylinder
)

