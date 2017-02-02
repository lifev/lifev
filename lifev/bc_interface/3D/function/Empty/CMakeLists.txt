SET(Empty_HEADERS
  Empty/BCInterfaceFunctionParserEmpty3D.hpp
  Empty/BCInterfaceFunctionUserDefinedEmpty3D.hpp
CACHE INTERNAL "")

SET(Empty_SOURCES
  Empty/BCInterfaceFunctionParserEmpty3D.cpp
  Empty/BCInterfaceFunctionUserDefinedEmpty3D.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${Empty_HEADERS})
  STRING(REPLACE "Empty/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/bc_interface/3D/function/Empty"
  COMPONENT ${PACKAGE_NAME}
)
