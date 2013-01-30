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
    @brief SolverOperator

    @author Umberto Villa <umberto.villa@gmail.com>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 03-09-2010
 */

#include <lifev/core/operator/SolverOperator.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"

#include <Teuchos_RCPBoostSharedPtrConversions.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wextra"

namespace LifeV
{

namespace Operators
{

SolverOperator::SolverOperator( boost::shared_ptr<Epetra_Comm> comm ):
	M_name( "SolverOperator" ),
	M_useTranspose( false ),
	M_lossOfAccuracy( undefined ),
    M_converged( undefined ),
    M_numIterations( 0 ),
    M_numCumulIterations( 0 ),
    M_tolerance( -1. ),
    M_printSubiterationCount( false ),
    M_comm( comm )
{ }

SolverOperator::~SolverOperator()
{
    M_prec.reset();
    M_oper.reset();
}

int SolverOperator::SetUseTranspose( bool useTranspose )
{
	M_useTranspose = useTranspose;

	int ierr( 0 );
	if( M_useTranspose )
		ierr = -1;

	return ierr;
}

void SolverOperator::setOperator( operatorPtr_Type _oper )
{
	ASSERT_PRE( _oper.get() != this, "Can't self assign" );
	ASSERT_PRE( _oper.get() != 0, "Can't assign a null pointer" );
	M_oper = _oper;
	doSetOperator();
}

void SolverOperator::setPreconditioner( operatorPtr_Type _prec )
{
	ASSERT_PRE( _prec.get() != this, "Self Assignment is forbidden" );
	ASSERT_PRE( _prec.get() != 0, "Can't assign a null pointer" );
	M_prec = _prec;
	doSetPreconditioner();
}

void SolverOperator::setParameters( const Teuchos::ParameterList& _pList )
{
	M_pList = Teuchos::rcp( new Teuchos::ParameterList( _pList ), true );
	doSetParameterList();
}

void SolverOperator::setTolerance( const Real& tolerance )
{
	M_tolerance = tolerance;
}

void SolverOperator::setUsedForPreconditioning( const bool& enable )
{
	M_printSubiterationCount = enable;
}

void SolverOperator::resetSolver()
{
    doResetSolver();
    M_prec.reset();
    M_oper.reset();
}

int SolverOperator::Apply( const vector_Type& X, vector_Type& Y ) const
{
	ASSERT_PRE( X.Map().SameAs( M_oper->OperatorDomainMap() ), "X and domain map do no coincide \n" );
	ASSERT_PRE( Y.Map().SameAs( M_oper->OperatorRangeMap() ) , "Y and range map do no coincide \n" );

	return M_oper->Apply( X, Y );
}

int SolverOperator::ApplyInverse( const vector_Type& X, vector_Type& Y ) const
{
	ASSERT_PRE( Y.Map().SameAs( M_oper->OperatorDomainMap() ), "Y and domain map do no coincide \n" );
	ASSERT_PRE( X.Map().SameAs( M_oper->OperatorRangeMap() ) , "X and range map do no coincide \n" );

	if ( M_useTranspose )
		return -1;
	int result = doApplyInverse( X, Y );

	M_numCumulIterations += M_numIterations;

    if( M_comm->MyPID() == 0 && M_printSubiterationCount )
    {
	    std::cout << "> " << numIterations() << " subiterations" << std::endl;
    }

	return result;
}

} // Namespace Operators

} // Namespace LifeV
