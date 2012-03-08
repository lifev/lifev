SET(solver_HEADERS
solver/OneDFSIData.hpp
solver/OneDFSIDefinitions.hpp
solver/OneDFSIFlux.hpp
solver/OneDFSIFluxLinear.hpp
solver/OneDFSIFluxNonLinear.hpp
solver/OneDFSIPhysics.hpp
solver/OneDFSIPhysicsLinear.hpp
solver/OneDFSIPhysicsNonLinear.hpp
solver/OneDFSISolver.hpp
solver/OneDFSISource.hpp
solver/OneDFSISourceLinear.hpp
solver/OneDFSISourceNonLinear.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
solver/OneDFSIData.cpp
solver/OneDFSIFluxLinear.cpp
solver/OneDFSIFluxNonLinear.cpp
solver/OneDFSIPhysicsLinear.cpp
solver/OneDFSIPhysicsNonLinear.cpp
solver/OneDFSISolver.cpp
solver/OneDFSISourceLinear.cpp
solver/OneDFSISourceNonLinear.cpp
CACHE INTERNAL "")
