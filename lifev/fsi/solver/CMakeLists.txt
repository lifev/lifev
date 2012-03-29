SET(solver_HEADERS
  solver/MonolithicRobinInterface.hpp
  solver/MonolithicBlock.hpp
  solver/MonolithicBlockComposedDNND.hpp
  solver/FSIMonolithicGI.hpp
  solver/FSIOperator.hpp
  solver/MonolithicBlockMatrixRN.hpp
  solver/FSISolver.hpp
  solver/MonolithicBlockComposedNN.hpp
  solver/MonolithicBlockComposedDND.hpp
  solver/FSIMonolithic.hpp
  solver/FSIExactJacobian.hpp
  solver/MonolithicBlockComposedDN.hpp
  solver/FSIData.hpp
  solver/MonolithicBlockComposed.hpp
  solver/HarmonicExtensionSolver.hpp
  solver/MonolithicBlockMatrix.hpp
  solver/FSIMonolithicGE.hpp
  solver/FSIFixedPoint.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/FSIMonolithicGI.cpp
  solver/MonolithicBlockComposedDND.cpp
  solver/MonolithicRobinInterface.cpp
  solver/MonolithicBlockMatrix.cpp
  solver/MonolithicBlockComposedNN.cpp
  solver/FSIMonolithic.cpp
  solver/MonolithicBlock.cpp
  solver/FSIOperator.cpp
  solver/FSIExactJacobian.cpp
  solver/MonolithicBlockComposed.cpp
  solver/MonolithicBlockComposedDNND.cpp
  solver/FSISolver.cpp
  solver/FSIMonolithicGE.cpp
  solver/FSIFixedPoint.cpp
  solver/MonolithicBlockComposedDN.cpp
  solver/FSIData.cpp
  solver/MonolithicBlockMatrixRN.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/fsi/solver"
  COMPONENT ${PACKAGE_NAME}
)
