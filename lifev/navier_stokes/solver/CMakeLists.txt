
INCLUDE(AddSubdirectories)
ADD_SUBDIRECTORIES(
  NavierStokesSolver
)

SET(solver_HEADERS
  solver/OseenAssembler.hpp
  solver/OseenData.hpp
  solver/StabilizationSD.hpp
  solver/StabilizationIP.hpp
  solver/OseenSolverShapeDerivative.hpp
  solver/OseenSolver.hpp
  ${NavierStokesSolver_HEADERS}
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/OseenData.cpp
  ${NavierStokesSolver_SOURCES}
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

FOREACH(INC ${NSLOCAL_HEADERS})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/navier_stokes/solver"
  COMPONENT ${PACKAGE_NAME}
)
