//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief This file contains the SolverBelosOperator class.

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 01-09-2011

    The class provides a way to use a the SolverBelos as an operator.
    This is particularly useful to apply "exactly the inverse of a
    block in a preconditioner.
 */

#include <SolverBelosOperator.hpp>
#include <life/lifecore/LifeDebug.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifefilters/GetPot.hpp>
#include <life/lifearray/VectorEpetra.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================

SolverBelosOperator::SolverBelosOperator( boost::shared_ptr<Epetra_Comm> comm ):
    Epetra_Operator(), M_comm(comm)
{

}

SolverBelosOperator::~SolverBelosOperator()
{
    M_solver.reset();
}

// ===================================================
// Methods
// ===================================================

Int
SolverBelosOperator::buildPreconditioner( operator_type& matrix, const list_Type& list )
{
    M_matrix = matrix;
    M_solver.reset( new linearSolver_type( M_comm ) );
    M_solver->setParameters( list );
    M_solver->setMatrix( M_matrix );

    return 0;
}

void
SolverBelosOperator::resetPreconditioner()
{
    M_solver.reset();
}

void
SolverBelosOperator::showMe( std::ostream& output ) const
{
    if( M_comm->MyPID() == 0) output << "SolverAmesosOperator:\n" << M_solver->parametersList() << std::endl;
}

// ===================================================
// Epetra Operator Interface Methods
// ===================================================
Int
SolverBelosOperator::SetUseTranspose( const bool /*useTranspose*/ )
{
    return 0;
}

bool
SolverBelosOperator::UseTranspose()
{
    return false;
}

Int
SolverBelosOperator::Apply( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    return M_matrix->matrixPtr()->Apply( X, Y );
}

Int
SolverBelosOperator::ApplyInverse( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    if( M_solver )
    {
        M_solver->setRightHandSide( X );
        M_solver->solve( Y );
    }
    return 0;
}

double
SolverBelosOperator::NormInf() const
{
    return 0;
}

const char*
SolverBelosOperator::Label() const
{
    return "SolverAmesosOperator";
}

bool
SolverBelosOperator::UseTranspose() const
{
    return false;
}

bool
SolverBelosOperator::HasNormInf() const
{
    return false;
}

const Epetra_Comm&
SolverBelosOperator::Comm() const
{
    return *M_comm;
}

const Epetra_Map&
SolverBelosOperator::OperatorRangeMap() const
{
    return M_matrix->matrixPtr()->OperatorRangeMap();
}

const Epetra_Map&
SolverBelosOperator::OperatorDomainMap() const
{
    return M_matrix->matrixPtr()->OperatorDomainMap();
}

} // namespace LifeV


