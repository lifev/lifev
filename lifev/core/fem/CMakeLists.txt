SET(fem_HEADERS
  fem/TimeAdvance.hpp
  fem/DOF.hpp
  fem/FEField.hpp
  fem/BCBase.hpp
  fem/PostProcessingBoundary.hpp
  fem/CurrentBoundaryFEBase.hpp
  fem/CurrentBoundaryFE.hpp
  fem/FEFunction.hpp
  fem/TimeData.hpp
  fem/TimeAdvanceData.hpp
  fem/GeometricMap.hpp
  fem/Assembly.hpp
  fem/ReferenceElement.hpp
  fem/FESpace.hpp
  fem/BCIdentifier.hpp
  fem/TimeAdvanceBDFVariableStep.hpp
  fem/BCManage.hpp
  fem/ReferenceFEHdiv.hpp
  fem/BCVector.hpp
  fem/DOFLocalPattern.hpp
  fem/QuadratureRule.hpp
  fem/BCManageNormal.hpp
  fem/CurrentFE.hpp
  fem/TimeAdvanceNewmark.hpp
  fem/ReferenceFEScalar.hpp
  fem/TimeAdvanceBDF.hpp
  fem/DOFInterface3Dto3D.hpp
  fem/QuadraturePoint.hpp
  fem/SobolevNorms.hpp
  fem/QuadratureRuleProvider.hpp
  fem/BCFunction.hpp
  fem/DOFInterface.hpp
  fem/ReferenceFEHybrid.hpp
  fem/DOFInterface3Dto2D.hpp
  fem/BCDataInterpolator.hpp
  fem/HyperbolicFluxNumerical.hpp
  fem/AssemblyElemental.hpp
  fem/ReferenceFE.hpp
  fem/BCHandler.hpp
CACHE INTERNAL "")

SET(fem_SOURCES
  fem/DOFInterface.cpp
  fem/CurrentBoundaryFEBase.cpp
  fem/DOFInterface3Dto3D.cpp
  fem/DOFInterface3Dto2D.cpp
  fem/BCFunction.cpp
  fem/AssemblyElemental.cpp
  fem/QuadratureRule.cpp
  fem/ReferenceFEScalar.cpp
  fem/BCVector.cpp
  fem/QuadraturePoint.cpp
  fem/BCIdentifier.cpp
  fem/CurrentBoundaryFE.cpp
  fem/BCBase.cpp
  fem/TimeData.cpp
  fem/TimeAdvanceData.cpp
  fem/BCManage.cpp
  fem/FEDefinitions.cpp
  fem/BCHandler.cpp
  fem/QuadratureRuleProvider.cpp
  fem/ReferenceFEHdiv.cpp
  fem/GeometricMap.cpp
  fem/DOFLocalPattern.cpp
  fem/ReferenceFEHybrid.cpp
  fem/CurrentFE.cpp
  fem/BCDataInterpolator.cpp
  fem/ReferenceElement.cpp
  fem/DOF.cpp
  fem/ReferenceFE.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${fem_HEADERS})
  STRING(REPLACE "fem/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/core/fem"
  COMPONENT ${PACKAGE_NAME}
)
