
INCLUDE(TribitsAddExecutable)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  SmoothAneurysm
  SOURCES main.cpp ud_functions.cpp flowConditions.cpp
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_exampleSmoothAneurysm
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)
