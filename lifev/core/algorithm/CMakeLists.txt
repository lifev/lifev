SET(HEADERS ${HEADERS}
ComposedOperator.hpp
EigenSolver.hpp
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
SolverAztecOO.hpp
)

SET(SOURCES ${SOURCES}
EigenSolver.cpp
PreconditionerAztecOO.cpp
PreconditionerComposed.cpp
Preconditioner.cpp
PreconditionerIfpack.cpp
PreconditionerML.cpp
SolverAmesos.cpp
SolverAztecOO.cpp
)
