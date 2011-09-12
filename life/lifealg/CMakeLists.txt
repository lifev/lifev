add_library (lifealg EigenSolver.cpp LinearEpetraOperatorBlock.cpp PreconditionerAztecOO.cpp PreconditionerComposed.cpp Preconditioner.cpp PreconditionerIfpack.cpp PreconditionerML.cpp SolverAmesos.cpp SolverAztecOO.cpp)
install (TARGETS lifealg DESTINATION lib)
install (FILES ComposedOperator.hpp EigenSolver.hpp LinearEpetraOperatorBlock.hpp LinearEpetraOperator.hpp NonLinearAitken.hpp NonLinearBrent.hpp NonLinearLineSearch.hpp NonLinearRichardson.hpp PreconditionerAztecOO.hpp PreconditionerComposed.hpp Preconditioner.hpp PreconditionerIfpack.hpp PreconditionerML.hpp SolverAmesos.hpp SolverAztecOO.hpp DESTINATION include/life/lifealg)

