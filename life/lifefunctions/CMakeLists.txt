add_library (lifefunctions BCInterfaceFunctionSolverDefined.cpp
OneDFSIFunctionSolverDefined.cpp RossEthierSteinmanDec.cpp
RossEthierSteinmanInc.cpp Womersley.cpp bessel/bessik.cpp
bessel/bessjy.cpp bessel/cbessik.cpp bessel/cbessjy.cpp bessel/gamma.cpp
)
install (TARGETS lifefunctions DESTINATION lib)
install (FILES BCInterfaceFactory.hpp BCInterfaceFunction.hpp BCInterfaceFunctionParserFile.hpp BCInterfaceFunctionParserFileSolver.hpp BCInterfaceFunctionParser.hpp BCInterfaceFunctionParserSolver.hpp BCInterfaceFunctionSolverDefined.hpp BCInterfaceFunctionUserDefined.hpp OneDFSIFunction.hpp OneDFSIFunctionSolverDefined.hpp RossEthierSteinmanDec.hpp RossEthierSteinmanInc.hpp Womersley.hpp DESTINATION include/life/lifefunctions)
install (FILES bessel.hpp DESTINATION include/life/lifefunctions/bessel)
