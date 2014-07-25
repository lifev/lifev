SET(solverIsotropic_HEADERS
  solver/isotropic/StructuralIsotropicConstitutiveLaw.hpp
  solver/isotropic/NeoHookeanMaterialNonLinear.hpp
  solver/isotropic/ExponentialMaterialNonLinear.hpp
  solver/isotropic/SecondOrderExponentialMaterialNonLinear.hpp
  solver/isotropic/VenantKirchhoffMaterialLinear.hpp
  solver/isotropic/VenantKirchhoffMaterialNonLinear.hpp
  solver/isotropic/VenantKirchhoffMaterialNonLinearPenalized.hpp
CACHE INTERNAL "")

SET(solverIsotropic_SOURCES
CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${solverIsotropic_HEADERS})
  STRING(REPLACE "solver/isotropic/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/structure/solver/isotropic"
  COMPONENT ${PACKAGE_NAME}_ISOTROPIC
)
