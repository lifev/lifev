INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  computePrincipalTensions
  SOURCES main.cpp ud_functions.cpp
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataPrincipalTensions
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CREATE_SYMLINK
)
