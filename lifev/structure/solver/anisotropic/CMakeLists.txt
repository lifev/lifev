SET(solverAnisotropic_HEADERS
  solver/anisotropic/StructuralAnisotropicConstitutiveLaw.hpp
  solver/anisotropic/HolzapfelMaterialNonLinear.hpp
  solver/anisotropic/HolzapfelGeneralizedMaterialNonLinear.hpp
  solver/anisotropic/AnisotropicMultimechanismMaterialNonLinear.hpp
  solver/anisotropic/DistributedHolzapfelMaterialNonLinear.hpp
CACHE INTERNAL "")

SET(solverAnisotropic_SOURCES
CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${solverAnisotropic_HEADERS})
  STRING(REPLACE "solver/anisotropic" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/structure/solver/anisotropic"
  COMPONENT ${PACKAGE_NAME}_ANISOTROPIC
)
