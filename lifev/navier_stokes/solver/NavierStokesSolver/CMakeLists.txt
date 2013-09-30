SET(NavierStokesSolver_HEADERS
    solver/NavierStokesSolver/AssemblyPolicyGeneralizedStokes.hpp
    solver/NavierStokesSolver/AssemblyPolicyStokes.hpp
    solver/NavierStokesSolver/AssemblyPolicyNavierStokesSemiImplicit.hpp
    solver/NavierStokesSolver/AssemblyPolicyNavierStokesNewton.hpp
    solver/NavierStokesSolver/AssemblyPolicyNavierStokesPicard.hpp
    solver/NavierStokesSolver/ExporterPolicyHDF5.hpp
    solver/NavierStokesSolver/ExporterPolicyNoExporter.hpp
    solver/NavierStokesSolver/InitPolicyInterpolation.hpp
    solver/NavierStokesSolver/InitPolicyProjection.hpp
    solver/NavierStokesSolver/InitPolicySolver.hpp
    solver/NavierStokesSolver/NavierStokesProblem.hpp
    solver/NavierStokesSolver/NavierStokesSolver.hpp
    solver/NavierStokesSolver/SolverPolicyLinearSolver.hpp
    solver/NavierStokesSolver/TimeIterationPolicyLinear.hpp
    solver/NavierStokesSolver/TimeIterationPolicyNonlinear.hpp
    solver/NavierStokesSolver/TimeIterationPolicyNonlinearIncremental.hpp
    CACHE INTERNAL "")

SET(NavierStokesSolver_SOURCES
    solver/NavierStokesSolver/AssemblyPolicyGeneralizedStokes.cpp
    solver/NavierStokesSolver/AssemblyPolicyStokes.cpp
    solver/NavierStokesSolver/AssemblyPolicyNavierStokesSemiImplicit.cpp
    solver/NavierStokesSolver/AssemblyPolicyNavierStokesNewton.cpp
    solver/NavierStokesSolver/AssemblyPolicyNavierStokesPicard.cpp
    solver/NavierStokesSolver/SolverPolicyLinearSolver.cpp
    CACHE INTERNAL "")

SET(NSLOCAL_HEADERS)
FOREACH(INC ${NavierStokesSolver_HEADERS})
  STRING(REPLACE "solver/NavierStokesSolver/" "" LOCAL_INC ${INC})
  SET(NSLOCAL_HEADERS ${NSLOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${NSLOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/navier_stokes/solver/NavierStokesSolver"
  COMPONENT ${PACKAGE_NAME}
)