SET(bc_HEADERS
  bc/BCInterface.hpp
  bc/BCInterfaceData.hpp
  bc/BCInterfaceDefinitions.hpp
CACHE INTERNAL "")

SET(bc_SOURCES
  bc/BCInterfaceData.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${bc_HEADERS})
  STRING(REPLACE "bc/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/bc_interface/core/bc"
  COMPONENT ${PACKAGE_NAME}
)
