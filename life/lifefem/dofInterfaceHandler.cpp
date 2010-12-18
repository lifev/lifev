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
    @brief Class for connecting the dof of a mesh (3D) and an interface (2D)
    that lives on the boundary of the mesh.

    @author V. Martin
    @date 00-02-2003

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#include <life/lifefem/dofInterfaceHandler.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

DofInterfaceHandler::DofInterfaceHandler( const UInt& NbNeigh ) :
        M_nbNeighbor( NbNeigh )
{
    M_neighborList.reserve( M_nbNeighbor );
    M_inputBCList.reserve( M_nbNeighbor );
    M_outputBCList.reserve( M_nbNeighbor );
    M_BCVectorList.reserve( M_nbNeighbor );

}

DofInterfaceHandler::DofInterfaceHandler( const UInt & NbNeigh, const LocalDofPattern& refFE, const Dof& dof1 ) :
        M_nbNeighbor( NbNeigh )
{
    M_neighborList.reserve( M_nbNeighbor );
    M_inputBCList.reserve( M_nbNeighbor );
    M_outputBCList.reserve( M_nbNeighbor );
    M_BCVectorList.reserve( M_nbNeighbor );

    for ( UInt i = 0 ; i < M_nbNeighbor ; ++i )
    {
        M_neighborList.push_back( dof_interface_type( new DofInterface3Dto2D( refFE, dof1 ) ) );
    }
}

// ===================================================
// Methods
// ===================================================

void DofInterfaceHandler::addNeighbor( const LocalDofPattern& refFE, const Dof& dof1 )
{
    ASSERT_PRE( M_nbNeighbor != M_neighborList.size(), "The list of neighbors in the Handler is full." );
    M_neighborList.push_back( dof_interface_type( new DofInterface3Dto2D( refFE, dof1 ) ) );
}


void DofInterfaceHandler::initVectors()
{
    ASSERT_PRE( M_nbNeighbor != M_inputBCList.size() || M_nbNeighbor == 0 ,
                "The list of InIBC Vectors in the Handler is full." );
    for ( UInt iter = 0 ; iter < M_nbNeighbor ; ++iter )
    {
        ASSERT_PRE( M_neighborList[ iter ]->finalized(),
                    "The DofInterface3Dto2D should be updated before calling initVectors (InIBC)." );

    }
    ASSERT_PRE( M_nbNeighbor != M_outputBCList.size() || M_nbNeighbor == 0 ,
                "The list of OutIBC Vectors in the Handler is full." );
    for ( UInt iter = 0 ; iter < M_nbNeighbor ; ++iter )
    {
        ASSERT_PRE( M_neighborList[ iter ]->finalized(),
                    "The DofInterface3Dto2D should be updated before calling initVectors (OutIBC)." );

        M_referenceToIndexMap[ M_neighborList[ iter ]->interfaceFlag() ] = iter;
    }
}


void DofInterfaceHandler::initBCVectorInterface()
{
    ASSERT_PRE( M_nbNeighbor != M_BCVectorList.size() || M_nbNeighbor == 0 ,
                "The list of BC Vectors in the Handler is full." );
    ASSERT_PRE( M_nbNeighbor == M_inputBCList.size(),
                "The list of InIBC Vectors in the Handler should be initialized before calling initBCVectorInterface." );
    ASSERT_PRE( M_nbNeighbor == M_outputBCList.size(),
                "The list of OutIBC Vectors in the Handler should be initialized before calling initBCVectorInterface." );
    for ( UInt iter = 0 ; iter < M_nbNeighbor ; ++iter )
    {
        ASSERT_PRE( M_neighborList[ iter ]->finalized(),
                    "The DofInterface3Dto2D should be updated before calling initBCVectorInterface." );

        BCVectorInterface::dofInterfacePtr_Type di = M_neighborList[ iter ];
        M_BCVectorList.push_back( BCVectorInterface( M_inputBCList[ iter ],
                                                     M_neighborList[ iter ]->nbInterfaceDof(), di ) );
    }

}

UInt DofInterfaceHandler::NbInterfaceUnknowns() const
{
    ASSERT_PRE( M_nbNeighbor == M_inputBCList.size(), "Some Vectors have not been added to the list (InIBC)." );
    ASSERT_PRE( M_nbNeighbor == M_outputBCList.size(), "Some Vectors have not been added to the list (OutIBC)." );
    UInt counter = 0;
    for ( UInt iter = 0 ; iter < M_nbNeighbor ; ++iter )
    {
        ASSERT_PRE( (int)M_neighborList[ iter ]->nbInterfaceDof() == M_inputBCList[ iter ].size() &&
                    (int)M_neighborList[ iter ]->nbInterfaceDof() == M_outputBCList[ iter ].size(),
                    "The IBC vectors must have the same size as the number of interface dof." );
        counter += M_neighborList[ iter ]->nbInterfaceDof();
    }
    return counter;
}

UInt DofInterfaceHandler::IndexOfInterfaceRef( const Int& interfref ) const
{
    std::map<Int, UInt>::const_iterator it = M_referenceToIndexMap.find( interfref );
    if ( it == M_referenceToIndexMap.end() )
    {
        std::ostringstream _err_msg;
        _err_msg << "Dof number " << interfref << " not found";
        ERROR_MSG( _err_msg.str().c_str() );
    }
    return it->second;
}

void DofInterfaceHandler::ClearSomeDofInterface3Dto2DList()
{
    for ( UInt iter = 0 ; iter < M_nbNeighbor ; ++iter )
    {
        ASSERT_PRE( M_neighborList[ iter ]->finalized(),
                    "The DofInterface3Dto2D should be updated before calling ClearDofInterface3Dto2DList()." );
        M_neighborList[ iter ]->clearLists();
    }
}

// ===================================================
// Operators
// ===================================================

//! Extracting neighbors from the list (starts from 0)
DofInterface3Dto2D& DofInterfaceHandler::operator[] ( const UInt& i )
{
    ASSERT_PRE( M_nbNeighbor == M_neighborList.size(), "Some neighbors have not been added to the list" );
    ASSERT_BD( i < M_nbNeighbor );
    return *M_neighborList[ i ];
}
const DofInterface3Dto2D& DofInterfaceHandler::operator[] ( const UInt& i ) const
{
    ASSERT_PRE( M_nbNeighbor == M_neighborList.size(), "Some neighbors have not been added to the list" );
    ASSERT_BD( i < M_nbNeighbor );
    return *M_neighborList[ i ];
}

// ===================================================
// Get Methods
// ===================================================

UInt DofInterfaceHandler::NbNeigh() const
{
    return M_neighborList.size();
}

const EpetraVector & DofInterfaceHandler::InIBC( const UInt & i ) const
{
    ASSERT_PRE( M_nbNeighbor == M_inputBCList.size(), "Some Vectors have not been added to the list (InIBC)." );
    ASSERT_BD( i < M_nbNeighbor );
    return M_inputBCList[ i ];
}
EpetraVector & DofInterfaceHandler::InIBC( const UInt & i )
{
    ASSERT_PRE( M_nbNeighbor == M_inputBCList.size(), "Some Vectors have not been added to the list (InIBC)." );
    ASSERT_BD( i < M_nbNeighbor );
    return M_inputBCList[ i ];
}

const EpetraVector & DofInterfaceHandler::InIBC_byRefInterf( const Int & refinterf ) const
{
    ASSERT_PRE( M_nbNeighbor == M_inputBCList.size(), "Some Vectors have not been added to the list (InIBC)." );
    UInt i = IndexOfInterfaceRef( refinterf );
    return M_inputBCList[ i ];
}
EpetraVector & DofInterfaceHandler::InIBC_byRefInterf( const Int & refinterf )
{
    ASSERT_PRE( M_nbNeighbor == M_inputBCList.size(), "Some Vectors have not been added to the list (InIBC)." );
    UInt i = IndexOfInterfaceRef( refinterf );
    return M_inputBCList[ i ];
}

const EpetraVector & DofInterfaceHandler::OutIBC( const UInt & i ) const
{
    ASSERT_PRE( M_nbNeighbor == M_outputBCList.size(), "Some Vectors have not been added to the list (OutIBC)." );
    ASSERT_BD( i < M_nbNeighbor );
    return M_outputBCList[ i ];
}
EpetraVector & DofInterfaceHandler::OutIBC( const UInt & i )
{
    ASSERT_PRE( M_nbNeighbor == M_outputBCList.size(), "Some Vectors have not been added to the list (OutIBC)." );
    ASSERT_BD( i < M_nbNeighbor );
    return M_outputBCList[ i ];
}

const EpetraVector & DofInterfaceHandler::OutIBC_byRefInterf( const Int & refinterf ) const
{
    ASSERT_PRE( M_nbNeighbor == M_outputBCList.size(), "Some Vectors have not been added to the list (OutIBC)." );
    UInt i = IndexOfInterfaceRef( refinterf );
    return M_outputBCList[ i ];
}
EpetraVector & DofInterfaceHandler::OutIBC_byRefInterf( const Int & refinterf )
{
    ASSERT_PRE( M_nbNeighbor == M_outputBCList.size(), "Some Vectors have not been added to the list (OutIBC)." );
    UInt i = IndexOfInterfaceRef( refinterf );
    return M_outputBCList[ i ];
}

const BCVectorInterface &DofInterfaceHandler:: BCvec( const UInt & i ) const
{
    ASSERT_PRE( M_nbNeighbor == M_BCVectorList.size(), "Some BC Vectors have not been added to the list" );
    return M_BCVectorList[ i ];
}
BCVectorInterface & DofInterfaceHandler::BCvec( const UInt & i )
{
    ASSERT_PRE( M_nbNeighbor == M_BCVectorList.size(), "Some BC Vectors have not been added to the list" );
    return M_BCVectorList[ i ];
}



}
