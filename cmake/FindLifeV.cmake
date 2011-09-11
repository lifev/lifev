#
# Find the UMFPACK includes and libraries
#

FIND_PATH(UMFPACK_INCLUDE_DIR umfpack.h amd.h UFconfig.h
  /usr/local/include
  /usr/include
  ${UMFPACK_ROOT}
  ${UMFPACK_ROOT}/include
)

find_library(UMFPACK_LIBRARY
  NAMES umfpack
  HINTS "/usr/local/lib"
        "/usr/lib"
        "${UMFPACK_ROOT}"
        "${UMFPACK_ROOT}/lib"
  )

find_library(AMD_LIBRARY
  NAMES amd
  HINTS "/usr/local/lib"
        "/usr/lib"
        "${UMFPACK_ROOT}"
        "${UMFPACK_ROOT}/lib"
  )

IF(UMFPACK_INCLUDE_DIR)
  IF(UMFPACK_LIBRARY)
    SET( UMFPACK_LIBRARIES ${UMFPACK_LIBRARY} ${AMD_LIBRARY})
    SET( UMFPACK_FOUND "YES" )
    MESSAGE ("-- Found UMFPACK")
  ENDIF(UMFPACK_LIBRARY)
ENDIF(UMFPACK_INCLUDE_DIR)
