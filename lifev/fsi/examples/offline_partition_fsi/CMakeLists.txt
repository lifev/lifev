INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  Offline_Partition_Write
  SOURCES main_write.cpp
  COMM mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_offline_partition_fsi
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(fluid_mesh_fsi_offline
  SOURCE_FILES tube20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
  CREATE_SYMLINK
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(solid_mesh_fsi_offline
  SOURCE_FILES vessel20.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/fsi/data/mesh/inria/
  CREATE_SYMLINK
)

TRIBITS_ADD_EXECUTABLE(
  Offline_Partition_Read
  SOURCES main_read.cpp
  COMM mpi
)
