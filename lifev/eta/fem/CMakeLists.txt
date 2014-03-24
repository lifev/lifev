SET(fem_HEADERS
  fem/ETCurrentFE.hpp
  fem/ETCurrentFE_FD3.hpp
  fem/ETCurrentFlag.hpp
  fem/ETCurrentBDFE.hpp
  fem/ETFESpace.hpp
  fem/MeshGeometricMap.hpp
  fem/QRAdapterBase.hpp
  fem/QRAdapterNeverAdapt.hpp
  fem/QuadratureBoundary.hpp
  fem/QuadratureRuleBoundary.hpp
  fem/LevelSetBDQRAdapter.hpp
CACHE INTERNAL "")

SET(fem_SOURCES
  fem/ETCurrentFE.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${fem_HEADERS})
  STRING(REPLACE "fem/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/eta/fem"
  COMPONENT ${PACKAGE_NAME}
)
