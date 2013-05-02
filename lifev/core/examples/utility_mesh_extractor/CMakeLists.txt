INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  UtilityMeshExtractor
  SOURCES main.cpp run_MeshExtractor.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_UtilityMeshExtractor
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_UtilityMeshExtractor
  SOURCE_FILES lumen.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/inria
)

