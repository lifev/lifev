add_library (lifearray MapEpetra.cpp MatrixElemental.cpp VectorElemental.cpp VectorEpetra.cpp)
install (TARGETS lifearray DESTINATION lib)
install (FILES ArraySimple.hpp MapEpetra.hpp MatrixContainer.hpp MatrixElemental.hpp MatrixEpetra.hpp MeshEntityContainer.hpp RNM.hpp RNMOperatorConstant.hpp RNMOperator.hpp RNMTemplate.hpp VectorContainer.hpp VectorElemental.hpp VectorEpetra.hpp VectorSmall.hpp DESTINATION include/life/lifearray)

