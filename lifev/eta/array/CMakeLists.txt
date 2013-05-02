SET(array_HEADERS
  array/ETMatrixElemental.hpp
  array/ETVectorElemental.hpp
  array/MatrixSmall.hpp
  array/OperationSmallAddition.hpp
  array/OperationSmallDivision.hpp
  array/OperationSmallDot.hpp
  array/OperationSmallProduct.hpp
  array/OperationSmallSubstraction.hpp
CACHE INTERNAL "")

SET(array_SOURCES ""
  array/ETMatrixElemental.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${array_HEADERS})
  STRING(REPLACE "array/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/eta/array"
  COMPONENT ${PACKAGE_NAME}
)
