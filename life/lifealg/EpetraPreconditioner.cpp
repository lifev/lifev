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
    @brief Epetra preconditioner

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 09-11-2006
 */

#include <lifeconfig.h>
#include "EpetraPreconditioner.hpp"

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
EpetraPreconditioner::EpetraPreconditioner( const comm_PtrType& comm ):
        M_precType              ( "EpetraPreconditioner" ),
        M_displayer             ( comm ),
        M_List                  (),
        M_preconditionerCreated ( false )
{
}

EpetraPreconditioner::EpetraPreconditioner( const EpetraPreconditioner& P, const comm_PtrType& comm ):
        M_precType              ( P.M_precType ),
        M_displayer             ( comm ),
        M_List                  ( P.getList() ),
        M_preconditionerCreated ( P.M_preconditionerCreated )
{
}

EpetraPreconditioner::~EpetraPreconditioner()
{
}

// ===================================================
// Methods
// ===================================================


// ===================================================
// Epetra Operator Interface Methods
// ===================================================
Int
EpetraPreconditioner::SetUseTranspose( const bool /*useTranspose=false*/ )
{
    assert( false );
    return 0;
}

Int
EpetraPreconditioner::Apply( const Epetra_MultiVector& /*X*/, Epetra_MultiVector& /*Y*/ ) const
{
    assert( false );
    return 0;
}

Int
EpetraPreconditioner::ApplyInverse( const Epetra_MultiVector& /*X*/, Epetra_MultiVector& /*Y*/ ) const
{
    assert( false );
    return 0;
}

bool
EpetraPreconditioner::UseTranspose()
{
    assert( false );
    return false;
}

const Epetra_Map&
EpetraPreconditioner::OperatorRangeMap() const
{
    assert( false );
    Epetra_Map *emptyMapPtr( NULL );
    return *emptyMapPtr;
}

const Epetra_Map&
EpetraPreconditioner::OperatorDomainMap() const
{
    assert( false );
    Epetra_Map *emptyMapPtr( NULL );
    return *emptyMapPtr;
}

// ===================================================
// Set Methods
// ===================================================
void
EpetraPreconditioner::setList( const list_Type& list )
{
    M_List = list;
}

void
EpetraPreconditioner::setSolver( SolverTrilinos& /*solver*/ )
{
    //assert( false );
}

// ===================================================
// Get Methods
// ===================================================
const bool&
EpetraPreconditioner::preconditionerCreated()
{
    return M_preconditionerCreated;
}

const EpetraPreconditioner::list_Type&
EpetraPreconditioner::getList() const
{
    return M_List;
}

EpetraPreconditioner::list_Type&
EpetraPreconditioner::list()
{
    return M_List;
}

} // namespace LifeV
