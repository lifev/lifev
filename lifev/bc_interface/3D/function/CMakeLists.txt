INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


SET(HEADERS "")
SET(SOURCES "")

ADD_SUBDIRECTORY(Empty)
SET(HEADERS ${HEADERS} ${Empty_HEADERS})
SET(SOURCES ${SOURCES} ${Empty_SOURCES})

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()
# This is the list of files excluded from the tarball
TRIBITS_EXCLUDE_FILES()

# Add library target
TRIBITS_ADD_LIBRARY(
  lifevbcinterface3Dempty
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS lifevbcinterface3D
)

IF (LifeV_ENABLE_NavierStokes)
SET(HEADERS "")
SET(SOURCES "")

ADD_SUBDIRECTORY(fluid)

SET(HEADERS ${HEADERS} ${fluid_HEADERS})
SET(SOURCES ${SOURCES} ${fluid_SOURCES})

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()
# This is the list of files excluded from the tarball
TRIBITS_EXCLUDE_FILES()

# Add library target
TRIBITS_ADD_LIBRARY(
  lifevbcinterface3Dfluid
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS lifevbcinterface3D
)
ENDIF()

IF (LifeV_ENABLE_Structure)
SET(HEADERS "")
SET(SOURCES "")

ADD_SUBDIRECTORY(solid)

SET(HEADERS ${HEADERS} ${solid_HEADERS})
SET(SOURCES ${SOURCES} ${solid_SOURCES})

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()
# This is the list of files excluded from the tarball
TRIBITS_EXCLUDE_FILES()

# Add library target
TRIBITS_ADD_LIBRARY(
  lifevbcinterface3Dsolid
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS lifevbcinterface3D
)
ENDIF()

IF (LifeV_ENABLE_FSI)
SET(HEADERS "")
SET(SOURCES "")

ADD_SUBDIRECTORY(fsi)

SET(HEADERS ${HEADERS} ${fsi_HEADERS})
SET(SOURCES ${SOURCES} ${fsi_SOURCES})

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()
# This is the list of files excluded from the tarball
TRIBITS_EXCLUDE_FILES()

# Add library target
TRIBITS_ADD_LIBRARY(
  lifevbcinterface3Dfsi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS lifevbcinterface3D
)
ENDIF()
