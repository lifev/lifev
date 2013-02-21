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
    @brief Implementation file for VectorBlockStructure

    @author Gwenol Grandperrin <gwenol.grandperrin@gmail.com>
    @date 21-08-2012
 */

#include <lifev/core/array/VectorBlockStructure.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================

VectorBlockStructure::
VectorBlockStructure()
    : M_blockSize(),
      M_blockFirstIndex(),
      M_totalSize( 0 )
{}

VectorBlockStructure::
VectorBlockStructure( const map_Type& map )
    : M_blockSize( 1, map.map( Unique )->NumGlobalElements() ),
      M_blockFirstIndex( 1, 0 ),
      M_totalSize( map.map( Unique )->NumGlobalElements() )
{}

VectorBlockStructure::
VectorBlockStructure( const mapVector_Type& mapVector )
    : M_blockSize( mapVector.nbMap() ),
      M_blockFirstIndex( mapVector.nbMap() )
{
    ASSERT( mapVector.nbMap() > 0 , "Map vector empty, impossible to construct a VectorBlockMonolithicEpetra!" );

	M_blockSize[0] = mapVector.mapSize( 0 );
	M_blockFirstIndex[0] = 0;

	M_totalSize = M_blockSize[0];

	for ( UInt i( 1 ); i < mapVector.nbMap(); ++i )
	{
		M_blockSize[i] = mapVector.mapSize( i );
		M_blockFirstIndex[i] = M_totalSize;

		M_totalSize += M_blockSize[i];
	}

}

VectorBlockStructure::
VectorBlockStructure( const VectorBlockStructure& blockStructure )
    : M_blockSize( blockStructure.M_blockSize ),
      M_blockFirstIndex( blockStructure.M_blockFirstIndex ),
      M_totalSize( blockStructure.M_totalSize )
{

}

// ===================================================
// Set Methods
// ===================================================

void
VectorBlockStructure::
setBlockStructure( const std::vector<UInt>& blockSizes )
{
    M_blockSize = blockSizes;

    M_blockFirstIndex.resize( M_blockSize.size() );

    M_totalSize = 0;

    for ( UInt i( 0 ); i< M_blockSize.size(); ++i )
    {
        M_blockFirstIndex[i] = M_totalSize;
        M_totalSize += M_blockSize[i];
    }
}

void
VectorBlockStructure::
setBlockStructure( const mapVector_Type& mapVector )
{
    ASSERT( mapVector.nbMap() > 0 , "Map vector empty, impossible to set the block structure" );

    M_blockSize.resize( mapVector.nbMap() );
    M_blockFirstIndex.resize( mapVector.nbMap() );

	M_totalSize = 0;

	for ( UInt i( 0 ); i < mapVector.nbMap(); ++i )
	{
		M_blockSize[i] = mapVector.mapSize( i );
		M_blockFirstIndex[i] = M_totalSize;
		M_totalSize += M_blockSize[i];
	}
}

void
VectorBlockStructure::
setBlockStructure( const VectorBlockStructure& blockStructure )
{
    UInt size( blockStructure.numBlocks() );
    M_blockSize.resize( size );
    M_blockFirstIndex.resize( size );

    M_totalSize = blockStructure.M_totalSize;

    for ( UInt i( 0 ); i < size; ++i )
    {
        M_blockSize[i] = blockStructure.M_blockSize[i];
        M_blockFirstIndex[i] = blockStructure.M_blockFirstIndex[i];
    }
}

} // Namespace LifeV
