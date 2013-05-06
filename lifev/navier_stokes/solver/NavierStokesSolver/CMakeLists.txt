# Define the header and sources files for the Navier-Stokes solver

set(NavierStokesSolver_SOURCES
    AssemblyPolicyGeneralizedStokes.cpp
    AssemblyPolicyStokes.cpp
    AssemblyPolicyNavierStokesSemiImplicit.cpp
    AssemblyPolicyNavierStokesNewton.cpp
    AssemblyPolicyNavierStokesPicard.cpp
    ExporterPolicyHDF5.cpp
    InitPolicyInterpolation.cpp
    NavierStokesProblem.cpp
    SolverPolicyLinearSolver.cpp
    CACHE INTERNAL ""
)
set(NavierStokesSolver_HEADERS
    AssemblyPolicyGeneralizedStokes.hpp
    AssemblyPolicyStokes.hpp
    AssemblyPolicyNavierStokesSemiImplicit.hpp
    AssemblyPolicyNavierStokesNewton.hpp
    AssemblyPolicyNavierStokesPicard.hpp
    ExporterPolicyHDF5.hpp
    ExporterPolicyNoExporter.hpp
    InitPolicyInterpolation.hpp
    InitPolicyProjection.hpp
    InitPolicySolver.hpp
    NavierStokesProblem.hpp
    NavierStokesSolver.hpp
    SolverPolicyLinearSolver.hpp
    TimeIterationPolicyLinear.hpp
    TimeIterationPolicyNonlinear.hpp
    TimeIterationPolicyNonlinearIncremental.hpp
    TimeIterationPolicyBlocksSolve.hpp
    CACHE INTERNAL ""
)

