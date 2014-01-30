SET(solid_HEADERS
  solid/BCInterfaceFunctionParserSolid3D.hpp
  solid/BCInterfaceFunctionParserSolverSolid3D.hpp
  solid/BCInterfaceFunctionSolverDefinedSolid3D.hpp
  solid/BCInterfaceFunctionUserDefinedSolid3D.hpp
CACHE INTERNAL "")

SET(solid_SOURCES
  solid/BCInterfaceFunctionParserSolid3D.cpp
  solid/BCInterfaceFunctionParserSolverSolid3D.cpp
  solid/BCInterfaceFunctionSolverDefinedSolid3D.cpp
  solid/BCInterfaceFunctionUserDefinedSolid3D.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solid_HEADERS})
  STRING(REPLACE "solid/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/bc_interface/3D/function/solid"
  COMPONENT ${PACKAGE_NAME}
)
