add_library ( lifearray MapEpetra.cpp
MatrixElemental.cpp
VectorEpetraStructuredView.cpp
VectorElemental.cpp
VectorEpetraStructured.cpp
VectorEpetra.cpp )
install ( TARGETS lifearray DESTINATION lib )
install ( FILES ArraySimple.hpp
MapEpetra.hpp
MapVector.hpp
MatrixContainer.hpp
MatrixElemental.hpp
MatrixEpetra.hpp
MatrixEpetraStructured.hpp
MatrixEpetraStructuredUtility.hpp
MatrixEpetraStructuredView.hpp
MeshEntityContainer.hpp
RNM.hpp
RNMOperatorConstant.hpp
RNMOperator.hpp
RNMTemplate.hpp
VectorContainer.hpp
VectorElemental.hpp
VectorEpetra.hpp
VectorEpetraStructured.hpp
VectorEpetraStructuredView.hpp
VectorSmall.hpp DESTINATION include/life/lifearray)
