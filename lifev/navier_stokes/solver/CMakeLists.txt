
SET(solver_HEADERS
  solver/OseenAssembler.hpp
  solver/OseenData.hpp
  solver/StabilizationSD.hpp
  solver/StabilizationIP.hpp
  solver/OseenSolverShapeDerivative.hpp
  solver/OseenSolver.hpp
  solver/StabilizationSUPG.hpp
  solver/StabilizationSUPGVMS.hpp
  solver/StabilizationVMSLES.hpp
  solver/NavierStokesSolver.hpp
  solver/SolverManager.hpp
  solver/aSIMPLE.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/OseenData.cpp
  solver/SolverManager.cpp
  solver/aSIMPLE.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/navier_stokes/solver"
  COMPONENT ${PACKAGE_NAME}
)

INCLUDE(AddSubdirectories)
ADD_SUBDIRECTORIES(
  NavierStokesSolver
)

SET(solver_HEADERS
  ${solver_HEADERS}
  ${NavierStokesSolver_HEADERS}
CACHE INTERNAL "")

SET(solver_SOURCES
  ${solver_SOURCES}
  ${NavierStokesSolver_SOURCES}
CACHE INTERNAL "")