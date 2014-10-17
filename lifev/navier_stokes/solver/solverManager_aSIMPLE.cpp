/*
 * OseenOperatorManagerLSC.cpp
 *
 *  Created on: Sep 30, 2010
 *      Author: uvilla
 */

#include <lifev/navier_stokes/solver/solverManager_aSIMPLE.hpp>
#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>
#include <lifev/core/linear_algebra/MLPreconditioner.hpp>
#include <lifev/core/linear_algebra/TwoLevelPreconditioner.hpp>

namespace LifeV
{

solverManager_aSIMPLE::solverManager_aSIMPLE():
        SolverManager(),
        M_aSIMPLEapproximatedMomentumOperator(new Operators::ApproximatedInvertibleRowMatrix),
        M_aSIMPLEapproximatedSchurComplementOperator(new Operators::ApproximatedInvertibleRowMatrix)
{
    M_approximatedMomentumOperator.reset(M_aSIMPLEapproximatedMomentumOperator);
    M_approximatedSchurComplementOperator.reset(M_aSIMPLEapproximatedSchurComplementOperator);
}

void solverManager_aSIMPLE::getMatrices(const matrixPtr_Type& F, const matrixPtr_Type& B, const matrixPtr_Type& Btranspose)
{
    M_F = F->matrixPtr();
    M_B = B->matrixPtr();
    M_Btranspose = Btranspose->matrixPtr();
}
    
void solverManager_aSIMPLE::updateApproximatedMomentumOperator( )
{
    M_displayer.leaderPrint( "\tupdate the approximate momentum operator...");
    LifeChrono chrono;
    chrono.start();
    M_aSIMPLEapproximatedMomentumOperator->SetRowMatrix(M_F);
    M_aSIMPLEapproximatedMomentumOperator->SetParameterList(*M_momentumList);
    M_aSIMPLEapproximatedMomentumOperator->Compute();
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );
}

void solverManager_aSIMPLE::updateApproximatedSchurComplementOperator( )
{
    M_displayer.leaderPrint( "\tupdate the approximate shur complement operator...");
    LifeChrono chrono;
    chrono.start();
    M_aSIMPLEapproximatedSchurComplementOperator->SetRowMatrix(M_F);
    M_aSIMPLEapproximatedSchurComplementOperator->SetParameterList(*M_schurComplementList);
    M_aSIMPLEapproximatedSchurComplementOperator->Compute();
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );
}

} /*end namespace */
