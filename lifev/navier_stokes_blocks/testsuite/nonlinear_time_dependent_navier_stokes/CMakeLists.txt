
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  challenge_steady
  SOURCES main.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_challenge_steady
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(solverOptions_challenge_steady
  CREATE_SYMLINK
  SOURCE_FILES solversOptionsFast.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

