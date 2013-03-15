
INCLUDE(TribitsAddExecutable)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  TensionsDisplCyl
  SOURCES main.cpp ud_functions.cpp
  COMM serial mpi
  )

# TRIBITS_COPY_FILES_TO_BINARY_DIR(dataCylindrical
#   SOURCE_FILES dataCylindrical
#   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
#   CREATE_SYMLINK
# )