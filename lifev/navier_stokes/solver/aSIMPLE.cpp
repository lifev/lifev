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
        SolverManager()
        //M_LSCapproximatedMomentumOperator(new Operators::ApproximatedInvertibleRowMatrix),
        //M_LSCapproximatedSchurComplementOperator(new Operators::LeastSquaresCommutator)
{
    //M_approximatedMomentumOperator.reset(M_LSCapproximatedMomentumOperator);
    //M_approximatedSchurComplementOperator.reset(M_LSCapproximatedSchurComplementOperator);
}

void aSIMPLE::updateApproximatedMomentumOperator( )
{
	/*
    M_LSCapproximatedMomentumOperator->SetRowMatrix(M_momentumMatrix);
    M_LSCapproximatedMomentumOperator->SetParameterList(*M_momentumOptions);
    M_LSCapproximatedMomentumOperator->Compute();
	*/
}

void aSIMPLE::updateApproximatedSchurComplementOperator( )
{
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
