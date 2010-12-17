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
    @brief EpetraMap

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 26-10-2006

    This class manages the distribution of elements of matrices or vectors on a parallel machine
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_Util.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/life.hpp>
#include <life/lifealg/EpetraMap.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
EpetraMap::EpetraMap():
    M_repeatedEpetraMap(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer(),
    M_commPtr()
{}

EpetraMap::EpetraMap( Int  numGlobalElements,
                      Int  numMyElements,
                      Int* myGlobalElements,
                      Int  indexBase,
                      const comm_ptrtype& commPtr ):
    M_repeatedEpetraMap(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer(),
    M_commPtr( commPtr )
{

    //Sort MyGlobalElements to avoid a bug in Trilinos (9?) when multiplying two matrices (A * B^T)
    std::sort ( myGlobalElements, myGlobalElements + numMyElements );

    createMap( numGlobalElements,
               numMyElements,
               myGlobalElements,
               indexBase,
               *commPtr );
}

EpetraMap::EpetraMap( const Int numGlobalElements,
                      const Int indexBase,
                      const comm_ptrtype& commPtr ) :
    M_repeatedEpetraMap(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer(),
    M_commPtr( commPtr )
{
    std::vector<Int> myGlobalElements( numGlobalElements );

    for ( Int i = 0; i < numGlobalElements; ++i )
        myGlobalElements[i] = i + indexBase;

    M_repeatedEpetraMap.reset( new Epetra_Map( -1, numGlobalElements, &myGlobalElements[0], indexBase, *commPtr ) );
    M_uniqueEpetraMap.reset( new Epetra_Map( numGlobalElements, indexBase, *commPtr ) );
}

EpetraMap::EpetraMap( const Int           size,
                      const comm_ptrtype& commPtr ):
    M_repeatedEpetraMap(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer(),
    M_commPtr( commPtr )
{
    Int numGlobalElements( size );
    Int numMyElements    ( numGlobalElements );
    std::vector<Int>  myGlobalElements( size );
    Int indexBase = 1;

    for ( Int i(0); i < numGlobalElements; ++i )
        myGlobalElements[i] = i + 1;

    M_repeatedEpetraMap.reset( new Epetra_Map( numGlobalElements,
                                               numMyElements,
                                               &myGlobalElements[0],
                                               indexBase,
                                               *commPtr ) );

    if ( commPtr->MyPID() != 0 ) numMyElements = 0;

    M_uniqueEpetraMap.reset( new Epetra_Map( numGlobalElements,
                                             numMyElements,
                                             &myGlobalElements[0],
                                             indexBase,
                                             *commPtr ) );
}

EpetraMap::EpetraMap( const map_type map ):
    M_repeatedEpetraMap( new map_type( map ) ),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer(),
    M_commPtr()
{
    uniqueMap();
}

EpetraMap::EpetraMap( const Epetra_BlockMap& blockMap, const Int offset, const Int maxId,
                      Int indexBase ) :
    M_repeatedEpetraMap(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer(),
    M_commPtr()
{

    if ( indexBase < 0 ) indexBase = blockMap.IndexBase();

    std::vector<Int> myGlobalElements;
    Int* sourceGlobalElements( blockMap.MyGlobalElements() );
    Int const startIdOrig( offset + indexBase );
    Int const endIdOrig  ( startIdOrig + maxId );
    const Int maxMyElements = std::min( maxId, blockMap.NumMyElements() );
    myGlobalElements.reserve( maxMyElements );

    //Sort MyGlobalElements to avoid a bug in Trilinos (9?) when multiplying two matrices (A * B^T)
    std::sort ( myGlobalElements.begin(), myGlobalElements.end() );

    // We consider that the source Map may not be ordered
    for ( Int i(0); i < blockMap.NumMyElements(); ++i )
        if ( sourceGlobalElements[i] < endIdOrig && sourceGlobalElements[i] >= startIdOrig )
            myGlobalElements.push_back( sourceGlobalElements[i] - offset );

    createMap( -1,
               myGlobalElements.size(),
               &myGlobalElements.front(),
               indexBase,
               blockMap.Comm() );
}

// ===================================================
// Operators
// ===================================================
EpetraMap &
EpetraMap::operator = ( const EpetraMap& epetraMap )
{

    if ( this != &epetraMap )
    {
        M_repeatedEpetraMap = epetraMap.M_repeatedEpetraMap;
        M_uniqueEpetraMap   = epetraMap.M_uniqueEpetraMap;
        M_exporter          = epetraMap.M_exporter;
        M_importer          = epetraMap.M_importer;
        M_commPtr           = epetraMap.M_commPtr;
    }

    return *this;
}


EpetraMap &
EpetraMap::operator += ( const EpetraMap& epetraMap )
{
    if ( ! epetraMap.getUniqueMap() )
        return *this;

    if ( ! this->getUniqueMap() )
    {
        this->operator = ( epetraMap );
        return *this;
    }

    Int*             pointer;
    std::vector<Int> map;

    pointer = getRepeatedMap()->MyGlobalElements();
    for ( Int ii = 0; ii < getRepeatedMap()->NumMyElements(); ++ii, ++pointer )
    {
        map.push_back( *pointer );
    }

    Int numGlobalElements = getUniqueMap()->NumGlobalElements()
                            + getUniqueMap()->IndexBase() - epetraMap.getUniqueMap()->IndexBase();

    pointer = epetraMap.getRepeatedMap()->MyGlobalElements();
    for (Int ii = 0; ii < epetraMap.getRepeatedMap()->NumMyElements(); ++ii, ++pointer)
    {
        map.push_back( *pointer + numGlobalElements );
    }

    Int IndexBase = getRepeatedMap()->IndexBase();

    M_repeatedEpetraMap.reset( new Epetra_Map(-1, map.size(), &map[0], IndexBase, epetraMap.getRepeatedMap()->Comm() ) );

    map.resize(0);
    pointer = getUniqueMap()->MyGlobalElements();

    for ( Int ii = 0; ii < getUniqueMap()->NumMyElements(); ++ii, ++pointer )
    {
        map.push_back( *pointer );
    }

    pointer = epetraMap.getUniqueMap()->MyGlobalElements();
    for ( Int ii = 0; ii < epetraMap.getUniqueMap()->NumMyElements(); ++ii, ++pointer )
    {
        map.push_back( *pointer + numGlobalElements );
    }

    M_uniqueEpetraMap.reset( new Epetra_Map( -1, map.size(), &map[0], IndexBase, epetraMap.getRepeatedMap()->Comm() ) );

    M_exporter.reset();
    M_importer.reset();

    return *this;
}

EpetraMap
EpetraMap::operator + ( const EpetraMap& epetraMap )
{
    EpetraMap map( *this );
    map += epetraMap;
    createImportExport();
    return map;
}

EpetraMap &
EpetraMap::operator += ( Int const size )
{
    EpetraMap  lagrMap( size, commPtr() );

    ASSERT( this->getUniqueMap(), "operator+=(const Int) works only for an existing EpetraMap" );

    this->operator+=( lagrMap );
    return *this;
}

EpetraMap
EpetraMap::operator +  ( Int const size )
{
    EpetraMap map( *this );
    map += size;
    createImportExport();
    return map;
}


// ===================================================
// Methods
// ===================================================
boost::shared_ptr<EpetraMap>
EpetraMap::createRootMap( Int const root )   const
{
    boost::shared_ptr<EpetraMap> rootMap( new EpetraMap( Epetra_Util::Create_Root_Map( *getUniqueMap(), root ) ) );
    return rootMap;
}

bool
EpetraMap::mapsAreSimilar( EpetraMap const& epetraMap ) const
{
    if ( this == &epetraMap )
        return true;

    return( getUniqueMap()->SameAs( *epetraMap.getUniqueMap() ) &&
            getRepeatedMap()->SameAs( *epetraMap.getRepeatedMap() ) );
}

void
EpetraMap::showMe( std::ostream& output ) const
{
    output << "showMe must be implemented for the EpetraMap class" << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
EpetraMap::map_ptrtype const &
EpetraMap::map( EpetraMapType mapType )   const
{
    switch ( mapType )
    {
    case Unique:
        return getUniqueMap();
    case Repeated:
        return getRepeatedMap();
    }
    return getUniqueMap();
}

Epetra_Export const&
EpetraMap::exporter()
{
    createImportExport();
    return **M_exporter;
}

Epetra_Import const&
EpetraMap::importer()
{
    createImportExport();
    return **M_importer;
}


// ===================================================
// Private Methods
// ===================================================
EpetraMap::EpetraMap( const EpetraMap& epetraMap ) :
    M_repeatedEpetraMap(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer()
{
    this->operator=( epetraMap );
}


void
EpetraMap::createMap( Int  numGlobalElements,
                      Int  numMyElements,
                      Int* myGlobalElements,
                      Int  indexBase,
                      const comm_type& comm )
{

    if ( numMyElements != 0 && myGlobalElements == 0 ) // linearMap
        M_repeatedEpetraMap.reset( new Epetra_Map( numGlobalElements,
                                                   numMyElements,
                                                   indexBase,
                                                   comm ) );
    else // classic LifeV map
        M_repeatedEpetraMap.reset( new Epetra_Map( numGlobalElements,
                                                   numMyElements,
                                                   myGlobalElements,
                                                   indexBase,
                                                   comm ) );

    uniqueMap();
}


void
EpetraMap::uniqueMap()
{
    M_uniqueEpetraMap.reset( new Epetra_Map( Epetra_Util::Create_OneToOne_Map ( *getRepeatedMap(), false ) ) );
    M_exporter.reset();
    M_importer.reset();
    return;
}

void
EpetraMap::createImportExport()
{

    if ( !getRepeatedMap() || !getUniqueMap() ) return;

    // The exporter is needed to import to a repeated vector
    if ( M_exporter.get() == 0 )
        M_exporter.reset( new boost::shared_ptr<Epetra_Export> );

    if ( M_exporter->get() == 0 )
        M_exporter->reset( new Epetra_Export( *getRepeatedMap(), *getUniqueMap() ) );

    if ( M_importer.get() == 0 )
        M_importer.reset( new boost::shared_ptr<Epetra_Import> );

    if ( M_importer->get() == 0 )
        M_importer->reset( new Epetra_Import( *getRepeatedMap(), *getUniqueMap() ) );

}

void
EpetraMap::bubbleSort(Epetra_IntSerialDenseVector& elements)
{
    Int hold;

    for ( Int pass(0); pass < elements.Length()-1; pass++ )
        for ( Int j(0); j < elements.Length()-1; j++ )
            if ( elements[j] > elements[j+1] )
            {
                hold          = elements[j];
                elements[j]   = elements[j+1];
                elements[j+1] = hold;
            }
}

void
EpetraMap::setUp( const RefFE&        refFE,
                  const comm_ptrtype& commPtr,
                  std::vector<Int>& repeatedNodeVector,
                  std::vector<Int>& repeatedEdgeVector,
                  std::vector<Int>& repeatedFaceVector,
                  std::vector<Int>& repeatedVolumeVector )
{
    Int indexBase = 1;

    if ( refFE.nbDofPerVertex() )
    {
        Int numNode = repeatedNodeVector.size();
        EpetraMap repeatedNodeMap( -1, numNode, &repeatedNodeVector[0], indexBase, commPtr );
        operator+=(repeatedNodeMap);
    }

    if ( refFE.nbDofPerEdge() )
    {
        Int numEdge = repeatedEdgeVector.size();
        EpetraMap repeatedEdgeMap( -1, numEdge, &repeatedEdgeVector[0], indexBase, commPtr );
        operator+=(repeatedEdgeMap);
    }

    if ( refFE.nbDofPerFace() )
    {
        Int numFace = repeatedFaceVector.size();
        EpetraMap repeatedFaceMap(-1, numFace, &repeatedFaceVector[0], indexBase, commPtr);
        operator+=( repeatedFaceMap );
    }

    if ( refFE.nbDofPerVolume() )
    {
        Int numElem = repeatedVolumeVector.size();
        EpetraMap repeatedElemMap( -1, numElem, &repeatedVolumeVector[0], indexBase, commPtr );
        operator+=( repeatedElemMap );
    }

    createImportExport();
}

} // end namespace LifeV

