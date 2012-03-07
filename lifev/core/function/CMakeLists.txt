INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
BCInterfaceFactory.hpp
BCInterfaceFunction.hpp
BCInterfaceFunctionParserFile.hpp
BCInterfaceFunctionParserFileSolver.hpp
BCInterfaceFunctionParser.hpp
BCInterfaceFunctionParserSolver.hpp
BCInterfaceFunctionSolverDefined.hpp
BCInterfaceFunctionUserDefined.hpp
)

SET(SOURCES ${SOURCES}
BCInterfaceFunctionSolverDefined.cpp
)
