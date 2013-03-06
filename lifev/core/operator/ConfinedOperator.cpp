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

ConfinedOperator::ConfinedOperator( boost::shared_ptr<Epetra_Comm> comm ):
        M_oper(),
        M_blockStructure(),
        M_blockIndex( 0 ),
        M_comm( comm ),
        M_map()
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
ConfinedOperator::setFullMap( const MapEpetra& map )
{
    M_map.reset( new Epetra_Map( *( map.map(Unique) ) ) );
}

void
ConfinedOperator::setBlockStructure( const blockStructure_Type& blockStructure )
{
    M_blockStructure.setBlockStructure( blockStructure );
}

void
ConfinedOperator::setBlockIndex( UInt index )
{
    ASSERT( M_blockStructure.numBlocks() > 0, "ConfinedOperator::setBlockIndex: Error: M_structure is not initialized" );
    ASSERT( index < M_blockStructure.numBlocks(), "ConfinedOperator::setBlockIndex: Error: index out of range" );
    M_blockIndex = index;
}

int
ConfinedOperator::Apply( const vector_Type& X, vector_Type& Y ) const
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::Apply: Error: M_oper pointer is null" );
    ASSERT( M_blockStructure.numBlocks() > 0, "ConfinedOperator::Apply: Error: M_structure is not initialized" );
    ASSERT( X.MyLength() == Y.MyLength(), "ConfinedOperator::Apply: Error: X and Y must have the same length" );

    int firstIndex = M_blockStructure.blockFirstIndex( M_blockIndex );
    Epetra_MultiVector xtmp( M_oper->OperatorDomainMap(), 1 );
    Epetra_MultiVector ytmp( M_oper->OperatorRangeMap() , 1 );

    // Extract the values from the vector
    const Int* gids         = M_oper->OperatorDomainMap().MyGlobalElements();
    const UInt numMyEntries = M_oper->OperatorDomainMap().NumMyElements();
    Int lid1;
    Int lid2;
    for ( UInt i = 0; i < numMyEntries; ++i )
    {
        lid1 = X.Map().LID( gids[i]+firstIndex );
        lid2 = M_oper->OperatorDomainMap().LID( gids[i] );
        ASSERT( ( lid2 >= 0 ) && ( lid1 >= 0 ), "ConfinedOperator::Apply: Error: lid < 0" );
        xtmp[0][lid2] = X[0][lid1];
    }

    // Apply the operator
    int result = M_oper->Apply( xtmp, ytmp );

    // Copy back the result in the Y vector;
    Y = X;
    const Int* gids2         = M_oper->OperatorRangeMap().MyGlobalElements();
    const UInt numMyEntries2 = M_oper->OperatorRangeMap().NumMyElements();
    for ( UInt i = 0; i < numMyEntries2; ++i )
    {
        lid1 = Y.Map().LID( gids2[i]+firstIndex );
        lid2 = M_oper->OperatorRangeMap().LID( gids2[i] );
        ASSERT( ( lid2 >= 0 ) && ( lid1 >= 0 ), "ConfinedOperator::Apply: Error: lid < 0" );
        Y[0][lid1] = ytmp[0][lid2];
    }

    return result;
}

int
ConfinedOperator::ApplyInverse( const vector_Type& X, vector_Type& Y ) const
{
    ASSERT( M_oper.get() != 0, "ConfinedOperator::ApplyInverse: Error: M_oper pointer is null" );
    ASSERT( M_blockStructure.numBlocks() > 0, "ConfinedOperator::ApplyInverse: Error: M_structure is not initialized" );
    ASSERT( X.MyLength() == Y.MyLength(), "ConfinedOperator::ApplyInverse: Error: X and Y must have the same length" );

    int firstIndex = M_blockStructure.blockFirstIndex( M_blockIndex );
    Epetra_MultiVector xtmp( M_oper->OperatorRangeMap() , 1 );
    Epetra_MultiVector ytmp( M_oper->OperatorDomainMap(), 1 );

    // Extract the values from the vector
    const Int* gids         = M_oper->OperatorRangeMap().MyGlobalElements();
    const UInt numMyEntries = M_oper->OperatorRangeMap().NumMyElements();
    Int lid1;
    Int lid2;
    for ( UInt i = 0; i < numMyEntries; ++i )
    {
        lid1 = X.Map().LID( gids[i]+firstIndex );
        lid2 = M_oper->OperatorRangeMap().LID( gids[i] );
        ASSERT( ( lid2 >= 0 ) && ( lid1 >= 0 ), "ConfinedOperator::ApplyInverse: Error: lid < 0" );
        xtmp[0][lid2] = X[0][lid1];
    }

    // Apply the operator
    int result = M_oper->ApplyInverse( xtmp, ytmp );

    // Copy back the result in the Y vector;
    Y = X;
    const Int* gids2         = M_oper->OperatorDomainMap().MyGlobalElements();
    const UInt numMyEntries2 = M_oper->OperatorDomainMap().NumMyElements();
    for ( UInt i = 0; i < numMyEntries2; ++i )
    {
        lid1 = Y.Map().LID( gids2[i]+firstIndex );
        lid2 = M_oper->OperatorDomainMap().LID( gids2[i] );
        ASSERT( ( lid2 >= 0 ) && ( lid1 >= 0 ), "ConfinedOperator::Apply: Error: lid < 0" );
        Y[0][lid1] = ytmp[0][lid2];
    }

    return result;
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
    ASSERT( M_blockStructure.numBlocks() > 0, "ConfinedOperator::OperatorDomainMap: Error: the structure is not known" );
    return *M_map;
}

const ConfinedOperator::map_Type&
ConfinedOperator::OperatorRangeMap() const
{
    ASSERT( M_blockStructure.numBlocks() > 0, "ConfinedOperator::OperatorRangeMap: Error: the structure is not known" );
    return *M_map;
}


} // Namespace Operators
} // Namespace LifeV

