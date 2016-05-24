
SET(solver_HEADERS
  solver/NavierStokesSolver.hpp
  solver/NavierStokesOperator.hpp
  solver/NavierStokesPreconditionerOperator.hpp
  solver/aSIMPLEOperator.hpp
  solver/aPCDOperator.hpp
  solver/Stabilization.hpp
#  solver/StabilizationSUPG.hpp
#  solver/StabilizationVMSLES.hpp
  solver/StabilizationSUPG_semi_implicit.hpp
  solver/StabilizationSUPG_semi_implicit_ale.hpp
  solver/StabilizationVMSLES_new.hpp
  solver/StabilizationVMSLES_semi_implicit.hpp
  solver/StabilizationSUPGALE.hpp
  solver/MultipleFluidsHandler.hpp
  solver/MultiFluidApplyOperator.hpp
  solver/DirichletNeumannPreconditionerMultiFluid.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/NavierStokesSolver.cpp
  solver/NavierStokesOperator.cpp
  solver/aSIMPLEOperator.cpp
  solver/aPCDOperator.cpp
#  solver/StabilizationSUPG.cpp
#  solver/StabilizationVMSLES.cpp
  solver/StabilizationSUPG_semi_implicit.cpp
  solver/StabilizationSUPG_semi_implicit_ale.cpp
  solver/StabilizationVMSLES_new.cpp
  solver/StabilizationVMSLES_semi_implicit.cpp
  solver/StabilizationSUPGALE.cpp
  solver/MultipleFluidsHandler.cpp
  solver/MultiFluidApplyOperator.cpp
  solver/DirichletNeumannPreconditionerMultiFluid.cpp
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