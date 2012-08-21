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
    @brief BelosOperator

    @author Gwenol Grandperrin <gwenol.grandperrin@gmail.com>

    @date 21-08-2012
 */

#include<lifev/core/operator/ConfinedOperator.hpp>

namespace LifeV
{
namespace Operators
{

ConfinedOperator::ConfinedOperator()
{

}

ConfinedOperator::~ConfinedOperator()
{

}


int
ConfinedOperator::SetUseTranspose( bool useTranspose )
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::SetUseTranspose: Error: M_oper pointer is null" );
    return M_oper->SetUseTranspose( useTranspose );
}

void
ConfinedOperator::setOperator( operatorPtr_Type oper )
{
    M_oper = oper;
}

void
ConfinedOperator::setStructure( blockStructurePtr_Type blockStructure )
{
    M_structure = blockStructure;
}

int
ConfinedOperator::Apply( const vector_Type& X, vector_Type& Y ) const
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::Apply: Error: M_oper pointer is null" );
    ASSERT( M_structure.get() != 0, "ConfinedOperator::Apply: Error: M_structure pointer is null" );
    return M_oper->ApplyInverse( X, Y );
}

int
ConfinedOperator::ApplyInverse( const vector_Type& X, vector_Type& Y ) const
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::ApplyInverse: Error: M_oper pointer is null" );
    ASSERT( M_structure.get() != 0, "ConfinedOperator::ApplyInverse: Error: M_structure pointer is null" );
    return M_oper->ApplyInverse( X, Y );
}

double
ConfinedOperator::NormInf() const
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::NormInf: Error: M_oper pointer is null" );
    return M_oper->NormInf();
}

const char*
ConfinedOperator::Label() const
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::Label: Error: M_oper pointer is null" );
    return M_oper->Label();
}

bool
ConfinedOperator::UseTranspose() const
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::UseTranspose: Error: M_oper pointer is null" );
    return M_oper->UseTranspose();
}

bool
ConfinedOperator::HasNormInf() const
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::HasNormInf: Error: M_oper pointer is null" );
    return M_oper->HasNormInf();
}

const ConfinedOperator::comm_Type&
ConfinedOperator::Comm() const
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::Comm: Error: M_oper pointer is null" );
    return M_oper->Comm();
}

const ConfinedOperator::map_Type&
ConfinedOperator::OperatorDomainMap() const
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::OperatorDomainMap: Error: M_oper pointer is null" );
    return M_oper->OperatorDomainMap();
}

const ConfinedOperator::map_Type&
ConfinedOperator::OperatorRangeMap() const
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::OperatorRangeMap: Error: M_oper pointer is null" );
    return M_oper->OperatorRangeMap();
}


} // Namespace Operators
} // Namespace LifeV

