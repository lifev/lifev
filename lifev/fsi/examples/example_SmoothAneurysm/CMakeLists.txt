
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  SmoothAneurysm
  SOURCES ud_functions.cpp flowConditions.cpp resistance.cpp main.cpp
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_exampleSmoothAneurysm
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)
