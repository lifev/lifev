SET(framework_HEADERS
  framework/MultiscaleCommunicatorsManager.hpp
  framework/MultiscaleDefinitions.hpp
  framework/MultiscaleGlobalData.hpp
  framework/MultiscaleInterface.hpp
  framework/MultiscaleSolver.hpp
CACHE INTERNAL "")

SET(framework_SOURCES
  framework/MultiscaleCommunicatorsManager.cpp
  framework/MultiscaleGlobalData.cpp
  framework/MultiscaleSolver.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${framework_HEADERS})
  STRING(REPLACE "framework/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/multiscale/framework"
  COMPONENT ${PACKAGE_NAME}
)
