SET(function_HEADERS
  function/BCInterfaceFunctionParser1D.hpp
  function/BCInterfaceFunctionParserSolver1D.hpp
  function/BCInterfaceFunctionSolverDefined1D.hpp
  function/BCInterfaceFunctionUserDefined1D.hpp
CACHE INTERNAL "")

SET(function_SOURCES
  function/BCInterfaceFunctionParser1D.cpp
  function/BCInterfaceFunctionParserSolver1D.cpp
  function/BCInterfaceFunctionSolverDefined1D.cpp
  function/BCInterfaceFunctionUserDefined1D.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${function_HEADERS})
  STRING(REPLACE "function/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/bc_interface/1D/function"
  COMPONENT ${PACKAGE_NAME}
)
