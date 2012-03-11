SET(solver_HEADERS
  solver/HeartMonodomainData.hpp
  solver/HeartStiffnessFibers.hpp
  solver/HeartIonicData.hpp
  solver/HeartBidomainData.hpp
  solver/HeartIonicSolver.hpp
  solver/HeartBidomainSolver.hpp
  solver/HeartFunctors.hpp
  solver/HeartMonodomainSolver.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/HeartIonicData.cpp
  solver/HeartMonodomainData.cpp
  solver/HeartBidomainData.cpp
  solver/HeartFunctors.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/heart/solver"
  COMPONENT ${PACKAGE_NAME}
)
