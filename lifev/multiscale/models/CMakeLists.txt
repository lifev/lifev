SET(models_HEADERS
  models/MultiscaleModel.hpp
  models/MultiscaleModelMultiscale.hpp
CACHE INTERNAL "")

SET(models_SOURCES
  models/MultiscaleModel.cpp
  models/MultiscaleModelMultiscale.cpp
CACHE INTERNAL "")





IF (LIFEV_HAS_ZERODIMENSIONAL)

SET(models0D_HEADERS
  models/MultiscaleModel0D.hpp
  models/MultiscaleModelWindkessel0D.hpp
CACHE INTERNAL "")
  
SET(models0D_SOURCES
  models/MultiscaleModel0D.cpp
  models/MultiscaleModelWindkessel0D.cpp
CACHE INTERNAL "")

ELSE()

SET(models0D_HEADERS "")
SET(models0D_SOURCES "")

ENDIF()





IF (LIFEV_HAS_ONEDFSI)

SET(models1D_HEADERS
  models/MultiscaleModelFSI1D.hpp
CACHE INTERNAL "")
  
SET(models1D_SOURCES
  models/MultiscaleModelFSI1D.cpp
CACHE INTERNAL "")
  
ELSE()

SET(models1D_HEADERS "")
SET(models1D_SOURCES "")

ENDIF()





IF (LIFEV_HAS_NAVIERSTOKES)

SET(modelsFluid3D_HEADERS
  models/MultiscaleModelFluid3D.hpp
CACHE INTERNAL "")
  
SET(modelsFluid3D_SOURCES
  models/MultiscaleModelFluid3D.cpp
CACHE INTERNAL "")
  
ELSE()

SET(modelsFluid3D_HEADERS "")
SET(modelsFluid3D_SOURCES "")
  
ENDIF()





IF (LIFEV_HAS_FSI)

SET(modelsFSI3D_HEADERS
  models/MultiscaleModelFSI3D.hpp
CACHE INTERNAL "")
  
SET(modelsFSI3D_SOURCES
  models/MultiscaleModelFSI3D.cpp
CACHE INTERNAL "")

ELSE()

SET(modelsFSI3D_HEADERS "")
SET(modelsFSI3D_SOURCES "")

ENDIF()





SET(models_HEADERS ${models_HEADERS} ${models0D_HEADERS})
SET(models_SOURCES ${models_SOURCES} ${models0D_SOURCES})

SET(models_HEADERS ${models_HEADERS} ${models1D_HEADERS})
SET(models_SOURCES ${models_SOURCES} ${models1D_SOURCES})

SET(models_HEADERS ${models_HEADERS} ${modelsFluid3D_HEADERS})
SET(models_SOURCES ${models_SOURCES} ${modelsFluid3D_SOURCES})

SET(models_HEADERS ${models_HEADERS} ${modelsFSI3D_HEADERS} CACHE INTERNAL "")
SET(models_SOURCES ${models_SOURCES} ${modelsFSI3D_SOURCES} CACHE INTERNAL "")






SET(LOCAL_HEADERS)
FOREACH(INC ${models_HEADERS})
  STRING(REPLACE "models/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/multiscale/models"
  COMPONENT ${PACKAGE_NAME}
)
