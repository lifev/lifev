SET(solver_HEADERS
  solver/ADRData.hpp
  solver/ADRAssemblerIP.hpp
  solver/ADRAssembler.hpp
  solver/HyperbolicSolver.hpp
  solver/HyperbolicData.hpp
CACHE INTERNAL "")


SET(solver_SOURCES
  solver/ADRData.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/core/solver"
  COMPONENT ${PACKAGE_NAME}
)
