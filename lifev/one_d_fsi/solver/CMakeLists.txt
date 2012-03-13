SET(solver_HEADERS
  solver/OneDFSIPhysicsNonLinear.hpp
  solver/OneDFSIFluxNonLinear.hpp
  solver/OneDFSISourceLinear.hpp
  solver/OneDFSIFlux.hpp
  solver/OneDFSIFluxLinear.hpp
  solver/OneDFSISource.hpp
  solver/OneDFSIPhysicsLinear.hpp
  solver/OneDFSISolver.hpp
  solver/OneDFSISourceNonLinear.hpp
  solver/OneDFSIDefinitions.hpp
  solver/OneDFSIData.hpp
  solver/OneDFSIPhysics.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/OneDFSISourceNonLinear.cpp
  solver/OneDFSIPhysicsLinear.cpp
  solver/OneDFSIData.cpp
  solver/OneDFSIPhysicsNonLinear.cpp
  solver/OneDFSISourceLinear.cpp
  solver/OneDFSIFluxLinear.cpp
  solver/OneDFSISolver.cpp
  solver/OneDFSIFluxNonLinear.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/one_d_fsi/solver"
  COMPONENT ${PACKAGE_NAME}
)
