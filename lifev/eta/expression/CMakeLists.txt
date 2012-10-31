SET(expression_HEADERS
	expression/EvaluationAddition.hpp
	expression/EvaluationDivI.hpp
	expression/EvaluationDivision.hpp
	expression/EvaluationDivJ.hpp
	expression/EvaluationDot.hpp
	expression/EvaluationDphiI.hpp
	expression/EvaluationDphiJ.hpp
	expression/EvaluationFunctor.hpp
	expression/EvaluationHK.hpp
	expression/EvaluationInterpolateGradient.hpp
	expression/EvaluationInterpolateValue.hpp
	expression/EvaluationMeas.hpp
	expression/EvaluationPhiI.hpp
	expression/EvaluationPhiJ.hpp
	expression/EvaluationPosition.hpp
	expression/EvaluationProduct.hpp
	expression/EvaluationScalar.hpp
	expression/EvaluationSubstraction.hpp
	expression/EvaluationVector.hpp
	expression/ExpressionAddition.hpp
	expression/ExpressionBase.hpp
	expression/ExpressionDivI.hpp
	expression/ExpressionDivision.hpp
	expression/ExpressionDivJ.hpp
	expression/ExpressionDot.hpp
	expression/ExpressionDphiI.hpp
	expression/ExpressionDphiJ.hpp
	expression/ExpressionFunctor.hpp
	expression/ExpressionHK.hpp
	expression/ExpressionInterpolateGradient.hpp
	expression/ExpressionInterpolateValue.hpp
	expression/ExpressionMeas.hpp
	expression/ExpressionPhiI.hpp
	expression/ExpressionPhiJ.hpp
	expression/ExpressionPosition.hpp
	expression/ExpressionProduct.hpp
	expression/ExpressionScalar.hpp
	expression/ExpressionSubstraction.hpp
	expression/ExpressionToEvaluation.hpp
	expression/ExpressionVector.hpp
	expression/Integrate.hpp
	expression/IntegrateMatrixElement.hpp
	expression/IntegrateValueElement.hpp
	expression/IntegrateVectorElement.hpp
	expression/RequestLoopElement.hpp
CACHE INTERNAL "")

SET(expression_SOURCES ""
	expression/ExpressionDivI.cpp
	expression/ExpressionDivJ.cpp
	expression/ExpressionDphiI.cpp
	expression/ExpressionDphiJ.cpp
	expression/ExpressionHK.cpp
	expression/ExpressionMeas.cpp
	expression/ExpressionPhiI.cpp
	expression/ExpressionPhiJ.cpp
	expression/ExpressionPosition.cpp
	expression/ExpressionScalar.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${expression_HEADERS})
  STRING(REPLACE "expression/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/eta/expression"
  COMPONENT ${PACKAGE_NAME}
)
