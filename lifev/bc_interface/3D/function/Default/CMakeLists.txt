SET(Default_HEADERS
  Default/BCInterfaceFunctionParserDefault3D.hpp
  Default/BCInterfaceFunctionUserDefinedDefault3D.hpp
CACHE INTERNAL "")

SET(Default_SOURCES
  Default/BCInterfaceFunctionParserDefault3D.cpp
  Default/BCInterfaceFunctionUserDefinedDefault3D.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${Default_HEADERS})
  STRING(REPLACE "Default/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/bc_interface/3D/function/Default"
  COMPONENT ${PACKAGE_NAME}
)
