SET(solver_HEADERS
  solver/MultiscaleAlgorithm.hpp
  solver/MultiscaleAlgorithmAitken.hpp
  solver/MultiscaleAlgorithmBroyden.hpp
  solver/MultiscaleAlgorithmExplicit.hpp
  solver/MultiscaleAlgorithmNewton.hpp
  solver/MultiscaleCommunicatorsManager.hpp
  solver/MultiscaleCoupling.hpp
  solver/MultiscaleCouplingBoundaryCondition.hpp
  solver/MultiscaleCouplingMeanNormalStress.hpp
  solver/MultiscaleCouplingMeanNormalStressArea.hpp
  solver/MultiscaleCouplingMeanNormalStressValve.hpp
  solver/MultiscaleCouplingMeanTotalNormalStress.hpp
  solver/MultiscaleCouplingMeanTotalNormalStressArea.hpp
  solver/MultiscaleGlobalData.hpp
  solver/MultiscaleDefinitions.hpp
  solver/MultiscaleInterface.hpp
  solver/MultiscaleModel.hpp
  solver/MultiscaleModel0D.hpp
  solver/MultiscaleModelFluid3D.hpp
  solver/MultiscaleModelFSI3D.hpp
  solver/MultiscaleModelFSI1D.hpp
  solver/MultiscaleModelWindkessel0D.hpp
  solver/MultiscaleModelMultiscale.hpp
  solver/MultiscaleSolver.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/MultiscaleAlgorithm.cpp
  solver/MultiscaleAlgorithmAitken.cpp
  solver/MultiscaleAlgorithmBroyden.cpp
  solver/MultiscaleAlgorithmExplicit.cpp
  solver/MultiscaleAlgorithmNewton.cpp
  solver/MultiscaleCommunicatorsManager.cpp
  solver/MultiscaleCoupling.cpp
  solver/MultiscaleCouplingBoundaryCondition.cpp
  solver/MultiscaleCouplingMeanNormalStress.cpp
  solver/MultiscaleCouplingMeanNormalStressArea.cpp
  solver/MultiscaleCouplingMeanNormalStressValve.cpp
  solver/MultiscaleCouplingMeanTotalNormalStress.cpp
  solver/MultiscaleCouplingMeanTotalNormalStressArea.cpp
  solver/MultiscaleGlobalData.cpp
  solver/MultiscaleModel.cpp
  solver/MultiscaleModel0D.cpp
  solver/MultiscaleModelFluid3D.cpp
  solver/MultiscaleModelFSI3D.cpp
  solver/MultiscaleModelFSI1D.cpp
  solver/MultiscaleModelWindkessel0D.cpp
  solver/MultiscaleModelMultiscale.cpp
  solver/MultiscaleSolver.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/multiscale/solver"
  COMPONENT ${PACKAGE_NAME}
)
