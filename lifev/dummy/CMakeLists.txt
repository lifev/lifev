
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)

#
# A) Define the package
#

TRIBITS_PACKAGE(Dummy)
SET(${PACKAGE_NAME}_VERSION "3.1")

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

SET(HEADERS "")
SET(SOURCES "")

ADD_SUBDIRECTORY(fem)
SET(HEADERS ${HEADERS} ${fem_HEADERS})
SET(SOURCES ${SOURCES} ${fem_SOURCES})

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()
# This is the list of files excluded from the tarball
TRIBITS_EXCLUDE_FILES()

# Add library target
TRIBITS_ADD_LIBRARY(
  lifevdummy
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
)

TRIBITS_ADD_TEST_DIRECTORIES(testsuite)
# TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

#
# Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
