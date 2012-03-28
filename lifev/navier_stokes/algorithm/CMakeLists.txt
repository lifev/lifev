SET(algorithm_HEADERS
  algorithm/PreconditionerLSC.hpp
  algorithm/PreconditionerPCD.hpp
  algorithm/PreconditionerSIMPLE.hpp
  algorithm/PreconditionerYosida.hpp
CACHE INTERNAL "")

SET(algorithm_SOURCES
  algorithm/PreconditionerLSC.cpp
  algorithm/PreconditionerPCD.cpp
  algorithm/PreconditionerSIMPLE.cpp
  algorithm/PreconditionerYosida.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${algorithm_HEADERS})
  STRING(REPLACE "algorithm/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/navier_stokes/algorithm"
  COMPONENT ${PACKAGE_NAME}
)
