SET(fluid_HEADERS
  fluid/BCInterfaceFunctionParserFluid3D.hpp
  fluid/BCInterfaceFunctionParserSolverFluid3D.hpp
  fluid/BCInterfaceFunctionUserDefinedFluid3D.hpp
CACHE INTERNAL "")

SET(fluid_SOURCES
  fluid/BCInterfaceFunctionParserFluid3D.cpp
  fluid/BCInterfaceFunctionParserSolverFluid3D.cpp
  fluid/BCInterfaceFunctionUserDefinedFluid3D.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${fluid_HEADERS})
  STRING(REPLACE "fluid/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/bc_interface/3D/function/fluid"
  COMPONENT ${PACKAGE_NAME}
)
