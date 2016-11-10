
SET(solver_HEADERS
  solver/NavierStokesSolverBlocks.hpp
  solver/NavierStokesOperator.hpp
  solver/NavierStokesPreconditionerOperator.hpp
  solver/aSIMPLEOperator.hpp
  solver/Stabilization.hpp
  solver/StabilizationSUPG.hpp
  solver/StabilizationSUPG_semi_implicit.hpp
  solver/StabilizationSUPG_semi_implicit_ale.hpp
  solver/StabilizationSUPGALE.hpp
  solver/StabilizationVMSLES_semi_implicit.hpp
  solver/FastAssemblerNS.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/NavierStokesSolverBlocks.cpp
  solver/NavierStokesOperator.cpp
  solver/aSIMPLEOperator.cpp
  solver/StabilizationSUPG.cpp
  solver/StabilizationSUPG_semi_implicit.cpp
  solver/StabilizationSUPG_semi_implicit_ale.cpp
  solver/StabilizationSUPGALE.cpp	
  solver/StabilizationVMSLES_semi_implicit.cpp
  solver/FastAssemblerNS.cpp
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