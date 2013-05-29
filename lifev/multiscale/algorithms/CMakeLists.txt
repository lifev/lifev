SET(algorithms_HEADERS
  algorithms/MultiscaleAlgorithm.hpp
  algorithms/MultiscaleAlgorithmAitken.hpp
  algorithms/MultiscaleAlgorithmExplicit.hpp
  algorithms/MultiscaleAlgorithmNewton.hpp
  algorithms/MultiscaleAlgorithmBroyden.hpp
CACHE INTERNAL "")

SET(algorithms_SOURCES
  algorithms/MultiscaleAlgorithm.cpp
  algorithms/MultiscaleAlgorithmAitken.cpp
  algorithms/MultiscaleAlgorithmExplicit.cpp
  algorithms/MultiscaleAlgorithmNewton.cpp
  algorithms/MultiscaleAlgorithmBroyden.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${algorithms_HEADERS})
  STRING(REPLACE "algorithms/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/multiscale/algorithms"
  COMPONENT ${PACKAGE_NAME}
)
