//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief AztecooOperator

    @author Umberto Villa <umberto.villa@gmail.com>

    @date 03-09-2010
 */

#include <lifev/core/operator/AztecooOperator.hpp>

namespace LifeV
{

namespace Operators
{

AztecooOperator::AztecooOperator():
		SolverOperator(),
		M_linSolver( new SolverType )
{
	M_name = "AztecOOOperator";
}

int
AztecooOperator::doApplyInverse( const vector_Type& X, vector_Type& Y ) const
{
    M_numIterations = 0;

	vector_Type Xcopy( X );
	Y.PutScalar( 0.0 );
	if( M_tolerance > 0 )
	    M_pList->sublist( "Trilinos: AztecOO List" ).set( "tol", M_tolerance );
	M_linSolver->SetParameters( M_pList->sublist( "Trilinos: AztecOO List" ) );
	M_linSolver->SetRHS( &Xcopy );
	M_linSolver->SetLHS( &Y );

    M_linSolver->SetUserOperator( M_oper.get() );

    if( M_prec.get() != 0 )
        M_linSolver->SetPrecOperator( (Epetra_Operator*) M_prec.get() );

	int maxIter( M_pList->sublist( "Trilinos: AztecOO List" ).get<int>( "max_iter" ) );
	double tol(  M_pList->sublist( "Trilinos: AztecOO List" ).get<double>( "tol" ) );

	// Solving the system
	int retValue = M_linSolver->Iterate(maxIter, tol);

    /* try to solve again (reason may be:
      -2 "Aztec status AZ_breakdown: numerical breakdown"
      -3 "Aztec status AZ_loss: loss of precision"
      -4 "Aztec status AZ_ill_cond: GMRES hessenberg ill-conditioned"
      This method was used in the old AztecOO solver.
    */
    if ( retValue <= -2 )
    {
        M_numIterations += M_linSolver->NumIters();
        retValue = M_linSolver->Iterate(maxIter, tol);
    }

	// Update the number of performed iterations
	M_numIterations += M_linSolver->NumIters();

	// Update of the status
	Real status[AZ_STATUS_SIZE];
	M_linSolver->GetAllAztecStatus( status );

	if ( status[AZ_why] == AZ_normal )
	{
		M_converged = yes;
	}
	else
	{
		M_converged = no;
	}

	if ( status[AZ_why] == AZ_loss )
	{
		M_lossOfAccuracy = yes;
	}
	else
	{
		M_lossOfAccuracy = no;
	}

	return retValue;
}

void
AztecooOperator::doSetOperator()
{

}

void
AztecooOperator::doSetPreconditioner()
{

}

void
AztecooOperator::doSetParameterList()
{
	M_linSolver->SetParameters( M_pList->sublist( "Trilinos: AztecOO List" ) );
}

void
AztecooOperator::doResetSolver()
{
    if( M_linSolver ) M_linSolver->DestroyPreconditioner();
}

} // Namespace Operators

} // Namespace LifeV
