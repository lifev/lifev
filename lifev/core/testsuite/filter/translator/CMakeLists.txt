INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE(
  ensightToVTK
  SOURCES ensightToVTK.cpp ../importExport/RossEthierSteinmanDec.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  vtkToEnsight
  SOURCES vtkToEnsight.cpp ../importExport/RossEthierSteinmanDec.cpp
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(dataTranslator
  SOURCE_FILES dataTranslator
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

