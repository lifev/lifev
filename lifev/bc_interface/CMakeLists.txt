
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)

#
# A) Define the package
#

TRIBITS_PACKAGE(BCInterface)
SET(${PACKAGE_NAME}_VERSION "3.1")

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()


ADD_SUBDIRECTORY(core)

IF (LifeV_ENABLE_ZeroDimensional)
ADD_SUBDIRECTORY(0D)
ENDIF ()

IF (LifeV_ENABLE_OneDFSI)
ADD_SUBDIRECTORY(1D)
ENDIF ()

ADD_SUBDIRECTORY(3D)

#
# Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
