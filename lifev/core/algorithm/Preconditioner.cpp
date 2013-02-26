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

#include <lifev/core/LifeV.hpp>
#include "Preconditioner.hpp"

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
Preconditioner::Preconditioner ( const commPtr_Type& comm ) :
    M_precType              ( "Preconditioner" ),
    M_displayer             ( comm ),
    M_list                  (),
    M_preconditionerCreated ( false )
{

}

Preconditioner::Preconditioner ( const Preconditioner& preconditioner, const commPtr_Type& comm ) :
    M_precType              ( preconditioner.M_precType ),
    M_displayer             ( comm ),
    M_list                  ( preconditioner.parametersList() ),
    M_preconditionerCreated ( preconditioner.M_preconditionerCreated )
{

}

Preconditioner::~Preconditioner()
{

}

// ===================================================
// Methods
// ===================================================


// ===================================================
// Epetra Operator Interface Methods
// ===================================================
Int
Preconditioner::SetUseTranspose ( const bool /*useTranspose=false*/ )
{
    assert ( false );
    return 0;
}

Int
Preconditioner::Apply ( const Epetra_MultiVector& /*vector1*/, Epetra_MultiVector& /*vector2*/ ) const
{
    assert ( false );
    return 0;
}

Int
Preconditioner::ApplyInverse ( const Epetra_MultiVector& /*vector1*/, Epetra_MultiVector& /*vector2*/ ) const
{
    assert ( false );
    return 0;
}

void
Preconditioner::showMe ( std::ostream& /*output*/ ) const
{
    assert ( false );
}

// ===================================================
// Set Methods
// ===================================================
void
Preconditioner::setParametersList ( const list_Type& list )
{
    M_list = list;
}

void
Preconditioner::setSolver ( SolverAztecOO& /*solver*/ )
{
    //assert( false );
}

// ===================================================
// Get Methods
// ===================================================
const bool&
Preconditioner::preconditionerCreated()
{
    return M_preconditionerCreated;
}

const Preconditioner::list_Type&
Preconditioner::parametersList() const
{
    return M_list;
}

Preconditioner::list_Type&
Preconditioner::parametersList()
{
    return M_list;
}

bool
Preconditioner::UseTranspose()
{
    assert ( false );
    return false;
}

const Epetra_Map&
Preconditioner::OperatorRangeMap() const
{
    assert ( false );
    Epetra_Map* emptyMapPtr ( NULL );
    return *emptyMapPtr;
}

const Epetra_Map&
Preconditioner::OperatorDomainMap() const
{
    assert ( false );
    Epetra_Map* emptyMapPtr ( NULL );
    return *emptyMapPtr;
}

} // namespace LifeV
