SET(couplings_HEADERS
  couplings/MultiscaleCoupling.hpp
  couplings/MultiscaleCouplingBoundaryCondition.hpp
  couplings/MultiscaleCouplingMeanNormalStress.hpp
  couplings/MultiscaleCouplingMeanNormalStressValve.hpp
  couplings/MultiscaleCouplingMeanTotalNormalStress.hpp
CACHE INTERNAL "")

SET(couplings_SOURCES
  couplings/MultiscaleCoupling.cpp
  couplings/MultiscaleCouplingBoundaryCondition.cpp
  couplings/MultiscaleCouplingMeanNormalStress.cpp
  couplings/MultiscaleCouplingMeanNormalStressValve.cpp
  couplings/MultiscaleCouplingMeanTotalNormalStress.cpp
CACHE INTERNAL "")





SET(couplingsArea_HEADERS "")
SET(couplingsArea_SOURCES "")

IF (LIFEV_HAS_FSI)

IF (LIFEV_HAS_ONEDFSI)

SET(couplingsArea_HEADERS
  couplings/MultiscaleCouplingMeanNormalStressArea.hpp
  couplings/MultiscaleCouplingMeanTotalNormalStressArea.hpp
CACHE INTERNAL "")
  
SET(couplingsArea_SOURCES
  couplings/MultiscaleCouplingMeanNormalStressArea.cpp
  couplings/MultiscaleCouplingMeanTotalNormalStressArea.cpp
CACHE INTERNAL "")

ENDIF()

ENDIF()



SET(couplings_HEADERS ${couplings_HEADERS} ${couplingsArea_HEADERS} CACHE INTERNAL "")
SET(couplings_SOURCES ${couplings_SOURCES} ${couplingsArea_SOURCES} CACHE INTERNAL "")




SET(LOCAL_HEADERS)
FOREACH(INC ${couplings_HEADERS})
  STRING(REPLACE "couplings/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/multiscale/couplings"
  COMPONENT ${PACKAGE_NAME}
)
