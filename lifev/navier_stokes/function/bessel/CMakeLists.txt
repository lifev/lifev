SET(function_bessel_HEADERS
function/bessel/bessel.hpp
CACHE INTERNAL "")

SET(function_bessel_SOURCES
function/bessel/bessik.cpp
function/bessel/bessjy.cpp
function/bessel/cbessik.cpp
function/bessel/cbessjy.cpp
function/bessel/gamma.cpp
CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${function_bessel_HEADERS})
      STRING(REPLACE "function/bessel/" "" LOCAL_INC ${INC})
      SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/navier_stokes/function/bessel"
  COMPONENT ${PACKAGE_NAME}
)