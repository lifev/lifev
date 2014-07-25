
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  RossEthierSteinman
  SOURCES main.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_rossEthierSteinman
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_rossEthierSteinman2
  SOURCE_FILES n2.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/navier_stokes/data/mesh/inria
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_rossEthierSteinman4
  SOURCE_FILES n4.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/navier_stokes/data/mesh/inria
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_rossEthierSteinman8
  SOURCE_FILES n8.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/navier_stokes/data/mesh/inria
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_rossEthierSteinman16
  SOURCE_FILES n16.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/navier_stokes/data/mesh/inria
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_rossEthierSteinman32
  SOURCE_FILES n32.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/navier_stokes/data/mesh/inria
)