
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  BasicExample
  SOURCES main.cpp heart.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_basic_examples_data
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_basic_examples_mesh
  CREATE_SYMLINK
  SOURCE_FILES lid16.mesh
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../data/mesh/inria
)

