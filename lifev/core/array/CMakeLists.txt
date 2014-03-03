SET(array_HEADERS
  array/EnumMapEpetra.hpp
  array/VectorEpetra.hpp
  array/MapVector.hpp
  array/VectorSmall.hpp
  array/RNMTemplate.hpp
  array/VectorBlockStructure.hpp
  array/MatrixBlockStructure.hpp
  array/MatrixEpetraStructuredUtility.hpp
  array/MatrixEpetraStructuredView.hpp
  array/MapEpetra.hpp
  array/MapEpetraData.hpp
  array/RNMOperatorConstant.hpp
  array/RNM.hpp
  array/VectorContainer.hpp
  array/MatrixElemental.hpp
  array/MatrixEpetra.hpp
  array/VectorEpetraStructured.hpp
  array/MatrixEpetraStructured.hpp
  array/VectorEpetraStructuredView.hpp
  array/ArraySimple.hpp
  array/VectorElemental.hpp
  array/RNMOperator.hpp
  array/MatrixContainer.hpp
  array/GhostHandler.hpp
CACHE INTERNAL "")

SET(array_SOURCES
  array/MatrixBlockStructure.cpp
  array/VectorElemental.cpp
  array/MatrixElemental.cpp
  array/VectorBlockStructure.cpp
  array/VectorEpetraStructuredView.cpp
  array/VectorEpetra.cpp
  array/MapEpetra.cpp
  array/VectorEpetraStructured.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${array_HEADERS})
  STRING(REPLACE "array/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/core/array"
  COMPONENT ${PACKAGE_NAME}
)
