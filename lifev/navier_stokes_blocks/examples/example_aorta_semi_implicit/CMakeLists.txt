
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  aorta_semi_implicit
  SOURCES main.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_aorta_semi_implicit
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(solverOptions_aorta_semi_implicit
  CREATE_SYMLINK
  SOURCE_FILES solversOptionsFast.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(solverOptionsLaplacian_aorta_semi_implicit
  CREATE_SYMLINK
  SOURCE_FILES SolverParamListLaplacian.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)