SET(function_HEADERS
  function/BCInterfaceFunctionParser0D.hpp
  function/BCInterfaceFunctionParserSolver0D.hpp
  function/BCInterfaceFunctionUserDefined0D.hpp
CACHE INTERNAL "")

SET(function_SOURCES
  function/BCInterfaceFunctionParser0D.cpp
  function/BCInterfaceFunctionParserSolver0D.cpp
  function/BCInterfaceFunctionUserDefined0D.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${function_HEADERS})
  STRING(REPLACE "function/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/bc_interface/0D/function"
  COMPONENT ${PACKAGE_NAME}
)
