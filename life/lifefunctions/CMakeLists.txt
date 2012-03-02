add_library ( lifefunctions KimMoin.cpp
OneDFSIFunctionSolverDefined.cpp
RossEthierSteinmanDec.cpp
bessel/cbessik.cpp
bessel/cbessjy.cpp
bessel/bessjy.cpp
bessel/bessik.cpp
bessel/gamma.cpp
Womersley.cpp
BCInterfaceFunctionSolverDefined.cpp
RossEthierSteinmanInc.cpp )
install ( TARGETS lifefunctions DESTINATION lib )
install ( FILES BCInterfaceFactory.hpp
BCInterfaceFunction.hpp
BCInterfaceFunctionParserFile.hpp
BCInterfaceFunctionParserFileSolver.hpp
BCInterfaceFunctionParser.hpp
BCInterfaceFunctionParserSolver.hpp
BCInterfaceFunctionSolverDefined.hpp
BCInterfaceFunctionUserDefined.hpp
KimMoin.hpp
OneDFSIFunction.hpp
OneDFSIFunctionSolverDefined.hpp
RossEthierSteinmanDec.hpp
RossEthierSteinmanInc.hpp
Womersley.hpp DESTINATION include/life/lifefunctions)
install ( FILES bessel/bessel.hpp DESTINATION include/life/lifefunctions/bessel)
