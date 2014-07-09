INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_0DMinimalModel
  SOURCES main.cpp
  ARGS "-i ParametersTNNP.xml -c"
  NUM_MPI_PROCS 1
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_0DMinimalModel_ParametersMCells
  CREATE_SYMLINK
  SOURCE_FILES ParametersMCells.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_0DMinimalModel_ParametersEndo
  CREATE_SYMLINK
  SOURCE_FILES ParametersEndo.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_0DMinimalModel_ParametersEpi
  CREATE_SYMLINK
  SOURCE_FILES ParametersEpi.xml 
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_0DMinimalModel_ParametersTNNP
  CREATE_SYMLINK
  SOURCE_FILES ParametersTNNP.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_0DMinimalModel_ParametersPB
  CREATE_SYMLINK
  SOURCE_FILES ParametersPB.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)
 



