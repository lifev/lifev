add_library ( lifealg SolverAmesos.cpp
PreconditionerML.cpp
PreconditionerComposed.cpp
LinearEpetraOperatorBlock.cpp
EigenSolver.cpp
PreconditionerIfpack.cpp
PreconditionerAztecOO.cpp
SolverAztecOO.cpp
Preconditioner.cpp )
install ( TARGETS lifealg DESTINATION lib )
install ( FILES ComposedOperator.hpp
EigenSolver.hpp
LinearEpetraOperatorBlock.hpp
LinearEpetraOperator.hpp
NonLinearAitken.hpp
NonLinearBrent.hpp
NonLinearLineSearch.hpp
NonLinearRichardson.hpp
PreconditionerAztecOO.hpp
PreconditionerComposed.hpp
Preconditioner.hpp
PreconditionerIfpack.hpp
PreconditionerML.hpp
SolverAmesos.hpp
SolverAztecOO.hpp DESTINATION include/life/lifealg)
