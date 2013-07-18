
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORY(multiscale)

IF (LifeV_ENABLE_OneDFSI)
ADD_SUBDIRECTORY(onedmodel)
ENDIF ()