SET(algorithm_HEADERS
algorithm/ComposedOperator.hpp
algorithm/EigenSolver.hpp
algorithm/NonLinearAitken.hpp
algorithm/NonLinearBrent.hpp
algorithm/NonLinearLineSearch.hpp
algorithm/NonLinearRichardson.hpp
algorithm/PreconditionerAztecOO.hpp
algorithm/PreconditionerComposed.hpp
algorithm/Preconditioner.hpp
algorithm/PreconditionerIfpack.hpp
algorithm/PreconditionerML.hpp
algorithm/SolverAmesos.hpp
algorithm/SolverAztecOO.hpp
CACHE INTERNAL "")

SET(algorithm_SOURCES
algorithm/EigenSolver.cpp
algorithm/PreconditionerAztecOO.cpp
algorithm/PreconditionerComposed.cpp
algorithm/Preconditioner.cpp
algorithm/PreconditionerIfpack.cpp
algorithm/PreconditionerML.cpp
algorithm/SolverAmesos.cpp
algorithm/SolverAztecOO.cpp
CACHE INTERNAL "")
