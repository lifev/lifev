SET(solver_HEADERS
  solver/VenantKirchhoffViscoelasticData.hpp
  solver/StructuralConstitutiveLawData.hpp
  solver/VenantKirchhoffViscoelasticSolver.hpp
  solver/StructuralConstitutiveLaw.hpp
  solver/StructuralOperator.hpp
  solver/NeoHookeanMaterialNonLinear.hpp
  solver/ExponentialMaterialNonLinear.hpp
  solver/VenantKirchhoffMaterialLinear.hpp
  solver/VenantKirchhoffMaterialNonLinear.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/VenantKirchhoffViscoelasticData.cpp
  solver/StructuralConstitutiveLawData.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/structure/solver"
  COMPONENT ${PACKAGE_NAME}
)
