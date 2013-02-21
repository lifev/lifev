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
    @brief Implementation file for VectorEpetraStructured

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 01 Jun 2011
 */

#include <lifev/core/array/VectorEpetraStructured.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

VectorEpetraStructured::
VectorEpetraStructured ( const map_type& map, const mapType_type& mapType )
    : VectorEpetra ( map, mapType ),
      M_blockStructure ( map )
{}

VectorEpetraStructured::
VectorEpetraStructured ( const mapVector_type& mapVector, const mapType_type& mapType )
    : VectorEpetra ( mapType ),
      M_blockStructure ( mapVector )
{
    ASSERT ( mapVector.nbMap() > 0 , "Map vector empty, impossible to construct a VectorBlockMonolithicEpetra!" );

    map_type myMap ( mapVector.totalMap() );

    // Set the global map
    this->setMap ( myMap );
}

VectorEpetraStructured::
VectorEpetraStructured ( const VectorEpetraStructured& vector )
    : VectorEpetra ( vector ),
      M_blockStructure ( vector.M_blockStructure )
{}

VectorEpetraStructured::
VectorEpetraStructured ( const VectorEpetraStructured& vector, const mapType_type& mapType )
    : VectorEpetra ( vector, mapType ),
      M_blockStructure ( vector.M_blockStructure )
{}

VectorEpetraStructured::
VectorEpetraStructured ( const VectorEpetraStructured& vector, const mapType_type& mapType, const combine_type& combineMode )
    : VectorEpetra ( vector, mapType, combineMode ),
      M_blockStructure ( vector.M_blockStructure )
{}

// ===================================================
// Set Methods
// ===================================================

void
VectorEpetraStructured::
setBlockStructure ( const std::vector<UInt>& blockSizes )
{
    M_blockStructure.setBlockStructure ( blockSizes );

    ASSERT ( M_blockStructure.totalSize() == static_cast<UInt> ( this->size() ), " Incompatible block structure (global size does not match) " );
}

void
VectorEpetraStructured::
setBlockStructure ( const mapVector_type& mapVector )
{
    M_blockStructure.setBlockStructure ( mapVector );

    ASSERT ( M_blockStructure.totalSize() == static_cast<UInt> ( this->size() ), " Incompatible block structure (global size does not match) " );
}

void
VectorEpetraStructured::
setBlockStructure ( const VectorBlockStructure& blockStructure )
{
    M_blockStructure.setBlockStructure ( blockStructure );
}

// ===================================================
// Get Methods
// ===================================================

UInt
VectorEpetraStructured::
blockSize ( const UInt& index ) const
{
    return M_blockStructure.blockSize ( index );
}

void
VectorEpetraStructured::
blockView ( const UInt& index, block_type& blockView )
{
    blockView.setup ( M_blockStructure.blockFirstIndex ( index ), M_blockStructure.blockSize ( index ), this );
}

VectorEpetraStructured::block_ptrType
VectorEpetraStructured::
block ( const UInt& index )
{
    block_ptrType vectorBlockView ( new block_type );

    vectorBlockView->setup ( M_blockStructure.blockFirstIndex ( index ), M_blockStructure.blockSize ( index ), this );

    return vectorBlockView;
}

VectorBlockStructure
VectorEpetraStructured::
blockStructure() const
{
    return M_blockStructure;
}


} // Namespace LifeV
