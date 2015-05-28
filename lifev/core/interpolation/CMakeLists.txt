SET(interpolation_HEADERS
  interpolation/RBFInterpolation.hpp
  interpolation/RBFlocallyRescaledScalar.hpp
CACHE INTERNAL "")

SET(interpolation_SOURCES

CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${interpolation_HEADERS})
  STRING(REPLACE "interpolation/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/core/interpolation"
  COMPONENT ${PACKAGE_NAME}
)
