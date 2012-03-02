add_library ( lifealg LinearEpetraOperatorBlock.cpp
SolverAztecOO.cpp
EigenSolver.cpp
Preconditioner.cpp
PreconditionerIfpack.cpp
PreconditionerComposed.cpp
PreconditionerML.cpp
SolverAmesos.cpp
PreconditionerAztecOO.cpp )
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
