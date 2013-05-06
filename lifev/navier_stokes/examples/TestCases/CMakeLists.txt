
INCLUDE(TribitsAddExecutable)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  NavierStokesFramework
  SOURCES main.cpp NavierStokesCavity.cpp NavierStokesEthierSteinman.cpp MeshUtility.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(datagetpot_nsframework
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataxml_nsframework
  SOURCE_FILES data.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

