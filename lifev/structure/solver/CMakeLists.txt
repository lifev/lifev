SET(solver_HEADERS
  solver/StructuralConstitutiveLawData.hpp
  solver/VenantKirchhoffViscoelasticSolver.hpp
  solver/StructuralOperator.hpp
  solver/WallTensionEstimator.hpp
  solver/WallTensionEstimatorCylindricalCoordinates.hpp
  solver/WallTensionEstimatorData.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/VenantKirchhoffViscoelasticData.cpp
  solver/StructuralConstitutiveLawData.cpp
  solver/WallTensionEstimatorData.cpp
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
