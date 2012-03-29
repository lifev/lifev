SET(fem_HEADERS
  fem/BCInterfaceData0D.hpp
  fem/BCInterface3D.hpp
  fem/BCInterface.hpp
  fem/BCInterfaceData.hpp
  fem/BCInterface0D.hpp
  fem/BCInterfaceDefinitions.hpp
  fem/BCInterfaceData1D.hpp
  fem/BCInterfaceData3D.hpp
  fem/BCInterface1D.hpp
CACHE INTERNAL "")

SET(fem_SOURCES
  fem/BCInterfaceData3D.cpp
  fem/BCInterfaceData.cpp
  fem/BCInterfaceData1D.cpp
  fem/BCInterfaceData0D.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${fem_HEADERS})
  STRING(REPLACE "fem/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/bc_interface/fem"
  COMPONENT ${PACKAGE_NAME}
)
