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
    @brief Implementation file for VectorBlockEpetra

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 01 Jun 2011
 */

#include <VectorBlockEpetra.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================

VectorBlockEpetra::
VectorBlockEpetra( const map_type& map, const mapType_type& mapType)
    : VectorEpetra(map,mapType),
      M_blockSize(1,map.map(Unique)->NumMyElements()),
      M_blockFirstIndex(1,0)
{}

VectorBlockEpetra::
VectorBlockEpetra( const mapVector_type& mapVector, const mapType_type& mapType)
    : VectorEpetra(mapType),
      M_blockSize(mapVector.nbMap()),
      M_blockFirstIndex(mapVector.nbMap())
{
    ASSERT( mapVector.nbMap() > 0 , "Map vector empty, impossible to construct a VectorBlockEpetra!");

    map_type myMap(mapVector.map(0));

	M_blockSize[0]=mapVector.map(0).map(Unique)->NumMyElements();
	M_blockFirstIndex[0]=0;

	UInt totalSize(M_blockSize[0]);

	for (UInt i(1); i<mapVector.nbMap(); ++i)
	{
		myMap += mapVector.map(i);
		M_blockSize[i]=mapVector.map(i).map(Unique)->NumMyElements();
		M_blockFirstIndex[i]=totalSize;

		totalSize+= M_blockSize[i];
	}

    // Set the global map
    this->setMap(myMap);
}

VectorBlockEpetra::
VectorBlockEpetra( const VectorBlockEpetra& vector)
    : VectorEpetra(vector),
      M_blockSize(vector.M_blockSize),
      M_blockFirstIndex(vector.M_blockFirstIndex)
{}

VectorBlockEpetra::
VectorBlockEpetra( const VectorBlockEpetra& vector, const mapType_type& mapType)
    : VectorEpetra(vector,mapType),
      M_blockSize(vector.M_blockSize),
      M_blockFirstIndex(vector.M_blockFirstIndex)
{}

VectorBlockEpetra::
VectorBlockEpetra( const VectorBlockEpetra& vector, const mapType_type& mapType, const combine_type& combineMode)
    : VectorEpetra(vector,mapType,combineMode),
      M_blockSize(vector.M_blockSize),
      M_blockFirstIndex(vector.M_blockFirstIndex)
{}



// ===================================================
// Operators
// ===================================================

// ===================================================
// Methods
// ===================================================

// ===================================================
// Set Methods
// ===================================================

void
VectorBlockEpetra::
setBlockStructure( const std::vector<UInt>& blockSizes)
{
    M_blockSize = blockSizes;

    M_blockFirstIndex.resize(M_blockSize.size());

    UInt currentSize(0);
    for (UInt i(0); i< M_blockSize.size(); ++i)
    {
        M_blockFirstIndex[i] = currentSize;
        currentSize += M_blockSize[i];
    }
}

// ===================================================
// Get Methods
// ===================================================

void
VectorBlockEpetra::
vectorBlockView( const UInt& index, block_type& blockView ) const
{
    blockView.setup( M_blockFirstIndex[index], M_blockSize[index], *this);
}

VectorBlockEpetra::block_ptrType
VectorBlockEpetra::
block( const UInt& index) const
{
    block_ptrType mbv(new block_type);

    mbv->setup( M_blockFirstIndex[index], M_blockSize[index], *this);

    return mbv;
}



// ===================================================
// Private Methods
// ===================================================

} // Namespace LifeV
