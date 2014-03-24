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
    @brief Implementation file for VectorBlockMonolithicEpetra

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 01 Jun 2011
 */

#include "VectorBlockMonolithicEpetra.hpp"

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

VectorBlockMonolithicEpetra::
VectorBlockMonolithicEpetra ( const map_type& map, const mapType_type& mapType)
    : VectorEpetra (map, mapType),
      M_blockSize (1, map.map (Unique)->NumGlobalElements() ),
      M_blockFirstIndex (1, 0)
{}

VectorBlockMonolithicEpetra::
VectorBlockMonolithicEpetra ( const mapVector_type& mapVector, const mapType_type& mapType)
    : VectorEpetra (mapType),
      M_blockSize (mapVector.nbMap() ),
      M_blockFirstIndex (mapVector.nbMap() )
{
    ASSERT ( mapVector.nbMap() > 0 , "Map vector empty, impossible to construct a VectorBlockMonolithicEpetra!");

    map_type myMap (mapVector.map (0) );

    M_blockSize[0] = mapVector.mapSize (0);
    M_blockFirstIndex[0] = 0;

    UInt totalSize (M_blockSize[0]);

    for (UInt i (1); i < mapVector.nbMap(); ++i)
    {
        myMap += mapVector.map (i);
        M_blockSize[i] = mapVector.mapSize (i);
        M_blockFirstIndex[i] = totalSize;

        totalSize += M_blockSize[i];
    }

    // Set the global map
    this->setMap (myMap);
}

VectorBlockMonolithicEpetra::
VectorBlockMonolithicEpetra ( const VectorBlockMonolithicEpetra& vector)
    : VectorEpetra (vector),
      M_blockSize (vector.M_blockSize),
      M_blockFirstIndex (vector.M_blockFirstIndex)
{}

VectorBlockMonolithicEpetra::
VectorBlockMonolithicEpetra ( const VectorBlockMonolithicEpetra& vector, const mapType_type& mapType)
    : VectorEpetra (vector, mapType),
      M_blockSize (vector.M_blockSize),
      M_blockFirstIndex (vector.M_blockFirstIndex)
{}

VectorBlockMonolithicEpetra::
VectorBlockMonolithicEpetra ( const VectorBlockMonolithicEpetra& vector, const mapType_type& mapType, const combine_type& combineMode)
    : VectorEpetra (vector, mapType, combineMode),
      M_blockSize (vector.M_blockSize),
      M_blockFirstIndex (vector.M_blockFirstIndex)
{}

// ===================================================
// Set Methods
// ===================================================

void
VectorBlockMonolithicEpetra::
setBlockStructure ( const std::vector<UInt>& blockSizes)
{
    M_blockSize = blockSizes;

    M_blockFirstIndex.resize (M_blockSize.size() );

    UInt currentSize (0);
    for (UInt i (0); i < M_blockSize.size(); ++i)
    {
        M_blockFirstIndex[i] = currentSize;
        currentSize += M_blockSize[i];
    }
}

void
VectorBlockMonolithicEpetra::
setBlockStructure ( const mapVector_type& mapVector)
{
    ASSERT ( mapVector.nbMap() > 0 , "Map vector empty, impossible to set the block structure");

    M_blockSize.resize (mapVector.nbMap() );
    M_blockFirstIndex.resize (mapVector.nbMap() );

    UInt totalSize (0);

    for (UInt i (0); i < mapVector.nbMap(); ++i)
    {
        M_blockSize[i] = mapVector.mapSize (i);
        M_blockFirstIndex[i] = totalSize;

        totalSize += M_blockSize[i];
    }

    ASSERT ( totalSize == this->size(), " Incompatible block structure (global size does not match) ");
}

// ===================================================
// Get Methods
// ===================================================

void
VectorBlockMonolithicEpetra::
blockView ( const UInt& index, block_type& blockView )
{
    ASSERT ( index < M_blockFirstIndex.size(), "Invalid block index");
    ASSERT ( index < M_blockSize.size(), "Invalid block index");

    blockView.setup ( M_blockFirstIndex[index], M_blockSize[index], this);
}

VectorBlockMonolithicEpetra::block_ptrType
VectorBlockMonolithicEpetra::
block ( const UInt& index)
{
    ASSERT ( index < M_blockFirstIndex.size(), "Invalid block index");
    ASSERT ( index < M_blockSize.size(), "Invalid block index");

    block_ptrType mbv (new block_type);

    mbv->setup ( M_blockFirstIndex[index], M_blockSize[index], this);

    return mbv;
}


} // Namespace LifeV
