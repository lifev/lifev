SET(solver_HEADERS
  solver/ZeroDimensionalSolver.hpp
  solver/ZeroDimensionalRythmosModelInterface.hpp
  solver/ZeroDimensionalDefinitions.hpp
  solver/ZeroDimensionalRythmosSolverInterface.hpp
  solver/ZeroDimensionalCircuitData.hpp
  solver/ZeroDimensionalData.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/ZeroDimensionalSolver.cpp
  solver/ZeroDimensionalRythmosModelInterface.cpp
  solver/ZeroDimensionalCircuitData.cpp
  solver/ZeroDimensionalRythmosSolverInterface.cpp
  solver/ZeroDimensionalData.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/zero_dimensional/solver"
  COMPONENT ${PACKAGE_NAME}
)
