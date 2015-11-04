SET(fsi_HEADERS
  fsi/BCInterfaceFunctionParserFSI3D.hpp
  fsi/BCInterfaceFunctionParserSolverFSI3D.hpp
  fsi/BCInterfaceFunctionSolverDefinedFSI3D.hpp
  fsi/BCInterfaceFunctionUserDefinedFSI3D.hpp
CACHE INTERNAL "")

SET(fsi_SOURCES
  fsi/BCInterfaceFunctionParserFSI3D.cpp
  fsi/BCInterfaceFunctionParserSolverFSI3D.cpp
  fsi/BCInterfaceFunctionSolverDefinedFSI3D.cpp
  fsi/BCInterfaceFunctionUserDefinedFSI3D.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${fsi_HEADERS})
  STRING(REPLACE "fsi/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/bc_interface/3D/function/fsi"
  COMPONENT ${PACKAGE_NAME}
)
