SET(solver_HEADERS
  solver/DarcySolverLinear.hpp
  solver/DarcySolverTransient.hpp
  solver/DarcySolverNonLinear.hpp
  solver/DarcySolverTransientNonLinear.hpp
  solver/DarcyData.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/DarcyDummy.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/darcy/solver"
  COMPONENT ${PACKAGE_NAME}
)
