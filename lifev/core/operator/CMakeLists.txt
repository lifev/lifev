SET(operator_HEADERS
  operator/AztecooOperator.hpp
  operator/BelosOperator.hpp
  operator/ConfinedOperator.hpp
  operator/LinearOperator.hpp
  operator/SolverOperator.hpp
CACHE INTERNAL "")

SET(operator_SOURCES
  operator/AztecooOperator.cpp
  operator/BelosOperator.cpp
  operator/ConfinedOperator.cpp
  operator/SolverOperator.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${operator_HEADERS})
  STRING(REPLACE "operator/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/core/operator"
  COMPONENT ${PACKAGE_NAME}
)
