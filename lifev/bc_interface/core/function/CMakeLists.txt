SET(function_HEADERS
  function/BCInterfaceFactory.hpp
  function/BCInterfaceFunction.hpp
  function/BCInterfaceFunctionParser.hpp
  function/BCInterfaceFunctionParserFile.hpp
  function/BCInterfaceFunctionParserFileSolver.hpp
  function/BCInterfaceFunctionParserSolver.hpp
  function/BCInterfaceFunctionSolverDefined.hpp
  function/BCInterfaceFunctionUserDefined.hpp
CACHE INTERNAL "")

SET(function_SOURCES
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${function_HEADERS})
  STRING(REPLACE "function/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/bc_interface/core/function"
  COMPONENT ${PACKAGE_NAME}
)
