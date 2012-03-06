
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)

#
# A) Define the package
#

TRIBITS_PACKAGE(Core)

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(LifeV_${PACKAGE_NAME}_ENABLE_SPIRIT_PARSER
  ENABLE_SPIRIT_PARSER
  "Enable Boost Spirit Parser"
  OFF )

FOREACH(TPL_NAME in ${Trilinos_TPL_LIST})
  IF(${TPL_NAME} STREQUAL "HDF5")
      SET(HAVE_HDF5 TRUE)
  ENDIF()
ENDFOREACH()

IF(TPL_Boost_ENABLED)
  SET(HAVE_BOOST TRUE)
ENDIF()

IF(TOUPPER(${CMAKE_BUILD_TYPE}) STREQUAL "DEBUG")
  SET(HAVE_LIFEV_DEBUG TRUE)
ENDIF()

IF(TPL_QHull_ENABLED)
  SET(HAVE_QHULL TRUE)
ENDIF()

FOREACH(TRILINOS_PACKAGE_NAME in ${Trilinos_PACKAGE_LIST})
  IF(${TRILINOS_PACKAGE_NAME} STREQUAL "RYTHMOS")
      SET(HAVE_TRILINOS_RYTHMOS TRUE)
  ENDIF()
ENDFOREACH()

FOREACH(Tpack ${NOX_LIBRARIES})
  IF(${Tpack} STREQUAL "noxthyra")
    SET(HAVE_NOX_THYRA TRUE)
  ENDIF()
ENDFOREACH()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_USING_HDF5
  HAVE_EPETRAEXT_HDF5
  "Enable EpetraExt interface support for HDF5. This interface requires an already installed HDF5 library\; the include, the library path and the library name must be specified as well."
  ${${PACKAGE_NAME}_ENABLE_HDF5_DEFAULT} )

#Adding a define for when PyTrilinos is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_PYTRILINOS "${${PROJECT_NAME}_ENABLE_PyTrilinos}")

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  doc/UserGuide
  doc/inout
  doc/matlab.README
  example/MapColoring/sample_map
  example/MapColoring/sample_matrix
  example/inout/build
  example/model_evaluator/GLpApp/Parallel2DMeshGeneratorFormat.pdf
  example/model_evaluator/GLpApp/README
  example/model_evaluator/GLpApp/data/square/generate-serial-meshes-1-2.sh
  example/model_evaluator/GLpApp/data/square/generate-serial-meshes-1-2.sh.out
  example/model_evaluator/GLpApp/data/square/square.1.2
  example/model_evaluator/GLpApp/data/square/square.1.2.000
  example/model_evaluator/GLpApp/data/square/square.1.2.001
  example/model_evaluator/GLpApp/data/square/square.1.2.edge
  example/model_evaluator/GLpApp/data/square/square.1.2.ele
  example/model_evaluator/GLpApp/data/square/square.1.2.epart.2
  example/model_evaluator/GLpApp/data/square/square.1.2.node
  example/model_evaluator/GLpApp/data/square/square.1.2.npart.2
  example/model_evaluator/GLpApp/data/square/square.1.2.poly
  example/model_evaluator/GLpApp/data/square/square.1.edge
  example/model_evaluator/GLpApp/data/square/square.1.ele
  example/model_evaluator/GLpApp/data/square/square.1.node
  example/model_evaluator/GLpApp/data/square/square.1.poly
  example/model_evaluator/GLpApp/data/square/square.2.edge
  example/model_evaluator/GLpApp/data/square/square.2.ele
  example/model_evaluator/GLpApp/data/square/square.2.node
  example/model_evaluator/GLpApp/data/square/square.2.poly
  example/model_evaluator/GLpApp/data/square/square.poly
  example/model_evaluator/GLpApp/from-triangle-to-serial-input-mesh.pl
  src/btf/pothen/btf_notes
  src/btf/pothen/main.f
  src/distdir
  src/transform/EpetraExt_Dirichlet_.*
  src/transform/EpetraExt_StaticCondensation_LinearProblem..*
  src/transform/EpetraExt_SubCopy_CrsMatrix..*
  src/zoltan/EpetraExt_ZoltanMpi.*
  test/Copy
  test/Makefile.template
  test/Zoltan/Dummy
  test/inout/build
  test/testAll.*
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
