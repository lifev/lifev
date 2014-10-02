
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  Turek_benchmark_unsteady
  SOURCES main.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_Turek_benchmark_unsteady
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)