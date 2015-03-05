
SET(solver_HEADERS
  solver/NavierStokesSolver.hpp
  solver/NavierStokesOperator.hpp
  solver/NavierStokesPreconditionerOperator.hpp
  solver/aSIMPLEOperator.hpp
  solver/aPCDOperator.hpp
  solver/StabilizationSUPG.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/NavierStokesSolver.cpp
  solver/NavierStokesOperator.cpp
  solver/aSIMPLEOperator.cpp
  solver/aPCDOperator.cpp
  solver/StabilizationSUPG.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/navier_stokes_blocks/solver"
  COMPONENT ${PACKAGE_NAME}
)