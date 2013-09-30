INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BasicTest
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data
  SOURCE_FILES RythmosAztecOOParams.xml RythmosAztecOOParamsLowsf.xml RythmosBelosParams.xml circuitFile.dat data.dat inputFile.dat
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)
