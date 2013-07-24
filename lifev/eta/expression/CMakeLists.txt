SET(expression_HEADERS
	expression/EvaluationAddition.hpp
	expression/EvaluationDivI.hpp
	expression/EvaluationDivision.hpp
	expression/EvaluationDivJ.hpp
	expression/EvaluationDot.hpp
	expression/EvaluationDphiI.hpp
	expression/EvaluationDphiJ.hpp
	expression/EvaluationEmult.hpp
	expression/EvaluationExtract1.hpp
	expression/EvaluationExtract2.hpp
	expression/EvaluationFunctor.hpp
	expression/EvaluationHK.hpp
	expression/EvaluationIfCrossed.hpp
	expression/EvaluationInterpolateGradient.hpp
	expression/EvaluationInterpolateValue.hpp
	expression/EvaluationMatrix.hpp
	expression/EvaluationMeas.hpp
	expression/EvaluationOuterProduct.hpp
	expression/EvaluationPhiI.hpp
	expression/EvaluationPhiJ.hpp
	expression/EvaluationPosition.hpp
	expression/EvaluationProduct.hpp
	expression/EvaluationNormal.hpp
	expression/EvaluationScalar.hpp
	expression/EvaluationSubstraction.hpp
	expression/EvaluationTranspose.hpp
	expression/EvaluationVector.hpp
	expression/ExpressionAddition.hpp
	expression/ExpressionBase.hpp
	expression/ExpressionDivI.hpp
	expression/ExpressionDivision.hpp
	expression/ExpressionDivJ.hpp
	expression/ExpressionDot.hpp
	expression/ExpressionDphiI.hpp
	expression/ExpressionDphiJ.hpp
	expression/ExpressionEmult.hpp
	expression/ExpressionExtract1.hpp
	expression/ExpressionExtract2.hpp
	expression/ExpressionFunctor.hpp
	expression/ExpressionHK.hpp
	expression/ExpressionIfCrossed.hpp
	expression/ExpressionInterpolateGradient.hpp
	expression/ExpressionInterpolateValue.hpp
	expression/ExpressionMatrix.hpp
	expression/ExpressionMeas.hpp
	expression/ExpressionOuterProduct.hpp
	expression/ExpressionPhiI.hpp
	expression/ExpressionPhiJ.hpp
	expression/ExpressionPosition.hpp
	expression/ExpressionNormal.hpp
	expression/ExpressionProduct.hpp
	expression/ExpressionScalar.hpp
	expression/ExpressionSubstraction.hpp	
	expression/ExpressionTranspose.hpp
	expression/ExpressionToEvaluation.hpp
	expression/ExpressionVector.hpp
	expression/Integrate.hpp
	expression/IntegrateMatrixElement.hpp
	expression/IntegrateMatrixElementLSAdapted.hpp
	expression/IntegrateMatrixFaceID.hpp
	expression/IntegrateMatrixFaceIDLSAdapted.hpp
	expression/IntegrateValueElement.hpp
	expression/IntegrateValueElementLSAdapted.hpp
	expression/IntegrateVectorElement.hpp
	expression/IntegrateVectorElementLSAdapted.hpp
	expression/IntegrateVectorFaceID.hpp
	expression/IntegrateVectorFaceIDLSAdapted.hpp
	expression/RequestLoopElement.hpp
	expression/RequestLoopFaceID.hpp
CACHE INTERNAL "")

SET(expression_SOURCES ""
	expression/EvaluationPhiI.cpp
	expression/EvaluationPhiJ.cpp
	expression/EvaluationScalar.cpp
	expression/ExpressionDivI.cpp
	expression/ExpressionDivJ.cpp
	expression/ExpressionDphiI.cpp
	expression/ExpressionDphiJ.cpp
	expression/ExpressionHK.cpp
	expression/ExpressionMeas.cpp
	expression/ExpressionPhiI.cpp
	expression/ExpressionPhiJ.cpp
	expression/ExpressionPosition.cpp
	expression/ExpressionNormal.cpp
	expression/ExpressionScalar.cpp
	expression/EvaluationPhiI.cpp
	expression/EvaluationPhiJ.cpp
	expression/EvaluationScalar.cpp
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
