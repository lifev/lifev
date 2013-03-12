SET(function_HEADERS
  function/BCInterfaceFunctionParser.hpp
  function/BCInterfaceFunctionSolverDefined.hpp
  function/BCInterfaceFunctionSolverDefinedFSI1D.hpp
  function/BCInterfaceFunctionSolverDefinedFSI3D.hpp
  function/BCInterfaceFunctionParserSolver.hpp
  function/BCInterfaceFunctionParserSolver0D.hpp
  function/BCInterfaceFunctionParserSolverFSI1D.hpp
  function/BCInterfaceFunctionParserSolverFSI3D.hpp
  function/BCInterfaceFunctionParserSolverFluid3D.hpp
  function/BCInterfaceFunctionUserDefined.hpp
  function/BCInterfaceFactory.hpp
  function/BCInterfaceFunctionParserFileSolver.hpp
  function/BCInterfaceFunction.hpp
  function/BCInterfaceFunctionParserFile.hpp
CACHE INTERNAL "")

SET(function_SOURCES
  function/BCInterfaceFunctionParserSolver0D.cpp
  function/BCInterfaceFunctionParserSolverFSI1D.cpp
  function/BCInterfaceFunctionParserSolverFSI3D.cpp
  function/BCInterfaceFunctionParserSolverFluid3D.cpp
  function/BCInterfaceFunctionSolverDefinedFSI1D.cpp
  function/BCInterfaceFunctionSolverDefinedFSI3D.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${function_HEADERS})
  STRING(REPLACE "function/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/bc_interface/function"
  COMPONENT ${PACKAGE_NAME}
)
