SET(algorithm_HEADERS
  algorithm/SolverAmesos.hpp
  algorithm/PreconditionerComposed.hpp
  algorithm/NonLinearLineSearch.hpp
  algorithm/PreconditionerAztecOO.hpp
  algorithm/NonLinearRichardson.hpp
  algorithm/ComposedOperator.hpp
  algorithm/SolverAztecOO.hpp
  algorithm/NonLinearAitken.hpp
  algorithm/PreconditionerIfpack.hpp
  algorithm/NonLinearBrent.hpp
  algorithm/Preconditioner.hpp
  algorithm/EigenSolver.hpp
  algorithm/LinearSolver.hpp
  algorithm/PreconditionerML.hpp
  algorithm/PreconditionerLinearSolver.hpp
CACHE INTERNAL "")

SET(algorithm_SOURCES
  algorithm/SolverAmesos.cpp
  algorithm/PreconditionerML.cpp
  algorithm/Preconditioner.cpp
  algorithm/PreconditionerAztecOO.cpp
  algorithm/PreconditionerComposed.cpp
  algorithm/PreconditionerIfpack.cpp
  algorithm/SolverAztecOO.cpp
  algorithm/EigenSolver.cpp
  algorithm/LinearSolver.cpp
  algorithm/PreconditionerLinearSolver.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${algorithm_HEADERS})
  STRING(REPLACE "algorithm/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/core/algorithm"
  COMPONENT ${PACKAGE_NAME}
)
