SET(linear_algebra_HEADERS
  linear_algebra/ApproximatedInvertibleRowMatrix.hpp
  linear_algebra/AztecooOperator.hpp
  linear_algebra/BelosOperator.hpp
  linear_algebra/BlockEpetra_Map.hpp
  linear_algebra/BlockEpetra_MultiVector.hpp
  linear_algebra/BlockOperator.hpp
  linear_algebra/CompositeOperator.hpp
  linear_algebra/IfpackPreconditioner.hpp
  linear_algebra/InvertibleOperator.hpp
  linear_algebra/LinearOperator.hpp
  linear_algebra/LumpedOperator.hpp
  linear_algebra/MLPreconditioner.hpp
  linear_algebra/RowMatrixPreconditioner.hpp
  linear_algebra/TwoLevelOperator.hpp
  linear_algebra/TwoLevelPreconditioner.hpp
CACHE INTERNAL "")

SET(linear_algebra_SOURCES
  linear_algebra/ApproximatedInvertibleRowMatrix.cpp
  linear_algebra/AztecooOperator.cpp
  linear_algebra/BelosOperator.cpp
  linear_algebra/BlockEpetra_Map.cpp
  linear_algebra/BlockEpetra_MultiVector.cpp
  linear_algebra/BlockOperator.cpp
  linear_algebra/CompositeOperator.cpp
  linear_algebra/IfpackPreconditioner.cpp
  linear_algebra/InvertibleOperator.cpp
  linear_algebra/LinearOperator.cpp
  linear_algebra/LumpedOperator.cpp
  linear_algebra/MLPreconditioner.cpp
  linear_algebra/TwoLevelOperator.cpp
  linear_algebra/TwoLevelPreconditioner.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${linear_algebra_HEADERS})
  STRING(REPLACE "linear_algebra/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/core/linear_algebra"
  COMPONENT ${PACKAGE_NAME}
)
