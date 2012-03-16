ADD_SUBDIRECTORY(bessel)

SET(function_HEADERS
  function/Womersley.hpp
  function/RossEthierSteinmanDec.hpp
  function/KimMoin.hpp
  function/RossEthierSteinmanInc.hpp
CACHE INTERNAL "")
  
SET(function_SOURCES
  function/Womersley.cpp
  function/RossEthierSteinmanInc.cpp
  function/KimMoin.cpp
  function/RossEthierSteinmanDec.cpp
CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${function_HEADERS})
      STRING(REPLACE "function/" "" LOCAL_INC ${INC})
      SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/navier_stokes/function"
  COMPONENT ${PACKAGE_NAME}
)

