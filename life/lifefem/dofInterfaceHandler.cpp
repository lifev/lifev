/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include "dofInterfaceHandler.hpp"

namespace LifeV
{
//! Default Constructor (call addNeigbor after for each interface)
/*! \param NbNeigh  number of neighbouring subdomains
*/
DofInterfaceHandler::DofInterfaceHandler( const UInt& NbNeigh ) :
        _nbNeigh( NbNeigh )
{
    _neighList.reserve( _nbNeigh );
    _InIBCList.reserve( _nbNeigh );
    _OutIBCList.reserve( _nbNeigh );
    _bcvList.reserve( _nbNeigh );

}

//! Constructor
/*! \param NbNeigh  number of neighbouring subdomains
    \param refFe the part of the reference FE that contains the dof patterns (nbDofPerEdge...)
    \param dof1 the Dof object of the mesh in which we want to make the computations
 */
DofInterfaceHandler::DofInterfaceHandler( const UInt & NbNeigh, const LocalDofPattern& refFE, const Dof& dof1 ) :
        _nbNeigh( NbNeigh )
{
    _neighList.reserve( _nbNeigh );
    _InIBCList.reserve( _nbNeigh );
    _OutIBCList.reserve( _nbNeigh );
    _bcvList.reserve( _nbNeigh );

    for ( UInt i = 0 ; i < _nbNeigh ; ++i )
    {
        _neighList.push_back( dof_interface_type( new DofInterface3Dto2D( refFE, dof1 ) ) );
    }
}

//! add a DofInterface3Dto2D to the list of neighbors
/*
    \param refFe the part of the reference FE that contains the dof patterns (nbDofPerEdge...)
    \param dof1 the Dof object of the mesh in which we want to make the computations
*/
void DofInterfaceHandler::addNeighbor( const LocalDofPattern& refFE, const Dof& dof1 )
{
    ASSERT_PRE( _nbNeigh != _neighList.size(), "The list of neighbors in the Handler is full." );
    _neighList.push_back( dof_interface_type( new DofInterface3Dto2D( refFE, dof1 ) ) );
}

//! creates the list of Vectors that store data
void DofInterfaceHandler::initVectors()
{
    ASSERT_PRE( _nbNeigh != _InIBCList.size() || _nbNeigh == 0 , "The list of InIBC Vectors in the Handler is full." );
    for ( UInt iter = 0 ; iter < _nbNeigh ; ++iter )
    {
        ASSERT_PRE( _neighList[ iter ]->finalized(), "The DofInterface3Dto2D should be updated before calling initVectors (InIBC)." );
        _InIBCList.push_back( Vector( _neighList[ iter ]->nbInterfaceDof() ) );

    }
    ASSERT_PRE( _nbNeigh != _OutIBCList.size() || _nbNeigh == 0 , "The list of OutIBC Vectors in the Handler is full." );
    for ( UInt iter = 0 ; iter < _nbNeigh ; ++iter )
    {
        ASSERT_PRE( _neighList[ iter ]->finalized(), "The DofInterface3Dto2D should be updated before calling initVectors (OutIBC)." );
        _OutIBCList.push_back( Vector( _neighList[ iter ]->nbInterfaceDof() ) );

        _indexInterfRefMap[ _neighList[ iter ]->InterfaceRef() ] = iter;
    }
}

//! creates the list of BC Vectors
void DofInterfaceHandler::initBCVectorInterface()
{
    ASSERT_PRE( _nbNeigh != _bcvList.size() || _nbNeigh == 0 , "The list of BC Vectors in the Handler is full." );
    ASSERT_PRE( _nbNeigh == _InIBCList.size(), "The list of InIBC Vectors in the Handler should be initialized before calling initBCVectorInterface." );
    ASSERT_PRE( _nbNeigh == _OutIBCList.size(), "The list of OutIBC Vectors in the Handler should be initialized before calling initBCVectorInterface." );
    for ( UInt iter = 0 ; iter < _nbNeigh ; ++iter )
    {
        ASSERT_PRE( _neighList[ iter ]->finalized(), "The DofInterface3Dto2D should be updated before calling initBCVectorInterface." );

        BCVectorInterface::dof_interface_type __di = _neighList[ iter ];
        _bcvList.push_back( BCVectorInterface( _InIBCList[ iter ], _neighList[ iter ]->nbInterfaceDof(), __di ) );
    }

}

//! How many neighbors stored?
UInt DofInterfaceHandler::NbNeigh() const
{
    //  std::cerr << "list size " << _neighList.size() << " and _nbNeigh " << _nbNeigh << std::endl;
    return _neighList.size();
}

//! Sum of the sizes of the Interface vectors
UInt DofInterfaceHandler::NbInterfaceUnknowns() const
{
    ASSERT_PRE( _nbNeigh == _InIBCList.size(), "Some Vectors have not been added to the list (InIBC)." );
    ASSERT_PRE( _nbNeigh == _OutIBCList.size(), "Some Vectors have not been added to the list (OutIBC)." );
    UInt counter = 0;
    for ( UInt iter = 0 ; iter < _nbNeigh ; ++iter )
    {
        ASSERT_PRE( _neighList[ iter ]->nbInterfaceDof() == _InIBCList[ iter ].size() &&
                    _neighList[ iter ]->nbInterfaceDof() == _OutIBCList[ iter ].size(),
                    "The IBC vectors must have the same size as the number of interface dof." );
        counter += _neighList[ iter ]->nbInterfaceDof();
    }
    return counter;
}

//! Extracting neighbors from the list (starts from 0)
DofInterface3Dto2D& DofInterfaceHandler::operator[] ( const UInt& i )
{
    ASSERT_PRE( _nbNeigh == _neighList.size(), "Some neighbors have not been added to the list" );
    ASSERT_BD( i >= 0 && i < _nbNeigh );
    return *_neighList[ i ];
}
const DofInterface3Dto2D& DofInterfaceHandler::operator[] ( const UInt& i ) const
{
    ASSERT_PRE( _nbNeigh == _neighList.size(), "Some neighbors have not been added to the list" );
    ASSERT_BD( i >= 0 && i < _nbNeigh );
    return *_neighList[ i ];
}


//! extracting a Vector in the _InIBCList list (starts from 0)
const Vector & DofInterfaceHandler::InIBC( const UInt & i ) const
{
    ASSERT_PRE( _nbNeigh == _InIBCList.size(), "Some Vectors have not been added to the list (InIBC)." );
    ASSERT_BD( i >= 0 && i < _nbNeigh );
    return _InIBCList[ i ];
}
Vector & DofInterfaceHandler::InIBC( const UInt & i )
{
    ASSERT_PRE( _nbNeigh == _InIBCList.size(), "Some Vectors have not been added to the list (InIBC)." );
    ASSERT_BD( i >= 0 && i < _nbNeigh );
    return _InIBCList[ i ];
}

//! extracting a Vector in the _InIBCList list (starts from 0)
//! using the reference of the interface
const Vector & DofInterfaceHandler::InIBC_byRefInterf( const Int & refinterf ) const
{
    ASSERT_PRE( _nbNeigh == _InIBCList.size(), "Some Vectors have not been added to the list (InIBC)." );
    UInt i = IndexOfInterfaceRef( refinterf );
    return _InIBCList[ i ];
}
Vector & DofInterfaceHandler::InIBC_byRefInterf( const Int & refinterf )
{
    ASSERT_PRE( _nbNeigh == _InIBCList.size(), "Some Vectors have not been added to the list (InIBC)." );
    UInt i = IndexOfInterfaceRef( refinterf );
    return _InIBCList[ i ];
}

//! extracting a Vector in the _OutIBCList list (starts from 0)
const Vector & DofInterfaceHandler::OutIBC( const UInt & i ) const
{
    ASSERT_PRE( _nbNeigh == _OutIBCList.size(), "Some Vectors have not been added to the list (OutIBC)." );
    ASSERT_BD( i >= 0 && i < _nbNeigh );
    return _OutIBCList[ i ];
}
Vector & DofInterfaceHandler::OutIBC( const UInt & i )
{
    ASSERT_PRE( _nbNeigh == _OutIBCList.size(), "Some Vectors have not been added to the list (OutIBC)." );
    ASSERT_BD( i >= 0 && i < _nbNeigh );
    return _OutIBCList[ i ];
}

//! extracting a Vector in the _OutIBCList list (starts from 0)
//! using the reference of the interface
const Vector & DofInterfaceHandler::OutIBC_byRefInterf( const Int & refinterf ) const
{
    ASSERT_PRE( _nbNeigh == _OutIBCList.size(), "Some Vectors have not been added to the list (OutIBC)." );
    UInt i = IndexOfInterfaceRef( refinterf );
    return _OutIBCList[ i ];
}
Vector & DofInterfaceHandler::OutIBC_byRefInterf( const Int & refinterf )
{
    ASSERT_PRE( _nbNeigh == _OutIBCList.size(), "Some Vectors have not been added to the list (OutIBC)." );
    UInt i = IndexOfInterfaceRef( refinterf );
    return _OutIBCList[ i ];
}


//! extracting a BCVector in the _bcvList list (starts from 0)
const BCVectorInterface &DofInterfaceHandler:: BCvec( const UInt & i ) const
{
    ASSERT_PRE( _nbNeigh == _bcvList.size(), "Some BC Vectors have not been added to the list" );
    return _bcvList[ i ];
}
BCVectorInterface & DofInterfaceHandler::BCvec( const UInt & i )
{
    ASSERT_PRE( _nbNeigh == _bcvList.size(), "Some BC Vectors have not been added to the list" );
    return _bcvList[ i ];
}


/*! This method returns the corresponding index number in the vectors
    (_neighList, _InIBCList...) living on the interfaces
    for a specific reference interface number.
  \param interfref : the reference of the interface.
*/
UInt DofInterfaceHandler::IndexOfInterfaceRef( const Int& interfref ) const
{
    std::map<Int, UInt>::const_iterator it = _indexInterfRefMap.find( interfref );
    if ( it == _indexInterfRefMap.end() )
        ERROR_MSG( "Dof number not found" );
    return it->second;
}

//! save some memory: destroy _neighList
//! use it safely!
//!(USE IT AFTER the construction of the interface mesh
//! and the _locDofMap in DofInterfaceBase)
void DofInterfaceHandler::ClearSomeDofInterface3Dto2DList()
{
    for ( UInt iter = 0 ; iter < _nbNeigh ; ++iter )
    {
        ASSERT_PRE( _neighList[ iter ]->finalized(), "The DofInterface3Dto2D should be updated before calling ClearDofInterface3Dto2DList()." );
        _neighList[ iter ]->ClearLists();
    }
}

}
