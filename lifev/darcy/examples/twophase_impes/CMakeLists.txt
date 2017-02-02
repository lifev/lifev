
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  TwoPhaseImpes
  SOURCES main.cpp impes.cpp user_fun.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_TwoPhaseImpes
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_TwoPhaseImpes_mesh
  CREATE_SYMLINK
  SOURCE_FILES mesh48603.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/darcy/data/mesh/inria
)
