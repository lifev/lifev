SET(solver_HEADERS
  solver/FSIHandler.hpp
  solver/FSIcouplingCE.hpp
  solver/ALESolver.hpp
  solver/FSIApplyOperator.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/FSIHandler.cpp
  solver/FSIcouplingCE.cpp
  solver/ALESolver.cpp
  solver/FSIApplyOperator.cpp	
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/fsi_blocks/solver"
  COMPONENT ${PACKAGE_NAME}
)
