/*
 * OseenOperatorManagerLSC.cpp
 *
 *  Created on: Sep 30, 2010
 *      Author: uvilla
 */

#include <lifev/navier_stokes/solver/aSIMPLE.hpp>
#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>
#include <lifev/core/linear_algebra/MLPreconditioner.hpp>
#include <lifev/core/linear_algebra/TwoLevelPreconditioner.hpp>

namespace LifeV
{

aSIMPLE::aSIMPLE():
        SolverManager(),
        M_aSIMPLEapproximatedMomentumOperator(new Operators::ApproximatedInvertibleRowMatrix),
        M_aSIMPLEapproximatedSchurComplementOperator(new Operators::ApproximatedInvertibleRowMatrix)
{
    M_approximatedMomentumOperator.reset(M_aSIMPLEapproximatedMomentumOperator);
    M_approximatedSchurComplementOperator.reset(M_aSIMPLEapproximatedSchurComplementOperator);
}

void aSIMPLE::getMatrices(const matrix_Type& F, const matrix_Type& B, const matrix_Type& Btranspose)
{
    M_F = F.matrixPtr();
    M_B = B.matrixPtr();
    M_Btranspose = Btranspose.matrixPtr();
}
    
void aSIMPLE::updateApproximatedMomentumOperator( )
{
    M_aSIMPLEapproximatedMomentumOperator->SetRowMatrix(M_F);
    M_aSIMPLEapproximatedMomentumOperator->SetParameterList(*M_momentumList);
    M_aSIMPLEapproximatedMomentumOperator->Compute();
}

void aSIMPLE::updateApproximatedSchurComplementOperator( )
{
    M_aSIMPLEapproximatedMomentumOperator->SetRowMatrix(M_F); // wrong
    M_aSIMPLEapproximatedMomentumOperator->SetParameterList(*M_schurComplementList);
    M_aSIMPLEapproximatedMomentumOperator->Compute();
    
	/*
    if( M_velocityLumpedMass.get() == 0 || M_recomputeConstantMatrices)
    {
        M_velocityLumpedMass.reset(new Operators::LumpedOperator);
        M_velocityLumpedMass->setUp(M_velocityMassMatrix);
        M_velocityLumpedMass->compute();
        M_recomputeConstantMatrices = true;
    }

    if( M_discreteLaplacian.get() == 0 || M_recomputeConstantMatrices)
    {
        Operators::DiscreteLaplacian * dl(new Operators::DiscreteLaplacian(M_displayer.comm() ) );
        dl->SetUp(M_pressureGradientMatrix, M_velocityLumpedMass);

        M_discreteLaplacian.reset(new Operators::ApproximatedInvertibleRowMatrix);
        M_discreteLaplacian->SetRowMatrix( dl->OperatorMatrix() );
        M_discreteLaplacian->SetParameterList( M_schurOptions->sublist("DiscreteLaplacian") );
        M_discreteLaplacian->Compute();
        delete dl;
    }

    M_LSCapproximatedSchurComplementOperator->SetUp(
            M_pressureGradientMatrix, M_velocityLumpedMass, M_discreteLaplacian);

    M_LSCapproximatedSchurComplementOperator->setC(M_momentumMatrix);
	*/
}

} /*end namespace */
