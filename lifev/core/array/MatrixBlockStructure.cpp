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
    @brief Implementation file for MatrixBlockStructure

    @author Gwenol Grandperrin <gwenol.grandperrin@gmail.com>
    @date 21-08-2012
 */

#include <lifev/core/array/MatrixBlockStructure.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

MatrixBlockStructure::MatrixBlockStructure()
    : M_rowsBlockStructure(), M_columnsBlockStructure()
{

}

MatrixBlockStructure::MatrixBlockStructure ( const map_Type& rowMap,
                                             const map_Type& columnMap )
    : M_rowsBlockStructure ( rowMap ), M_columnsBlockStructure ( columnMap )
{

}

MatrixBlockStructure::MatrixBlockStructure ( const map_Type& map )
    : M_rowsBlockStructure ( map ), M_columnsBlockStructure ( map )
{

}

MatrixBlockStructure::MatrixBlockStructure ( const mapVector_Type& rowMapVector,
                                             const mapVector_Type& columnMapVector )
    : M_rowsBlockStructure ( rowMapVector ), M_columnsBlockStructure (columnMapVector )
{

}

MatrixBlockStructure::MatrixBlockStructure ( const mapVector_Type& mapVector )
    : M_rowsBlockStructure ( mapVector ), M_columnsBlockStructure ( mapVector )
{

}

MatrixBlockStructure::MatrixBlockStructure ( const VectorBlockStructure& rowsBlockStructure,
                                             const VectorBlockStructure& columnsBlockStructure )
    : M_rowsBlockStructure ( rowsBlockStructure ), M_columnsBlockStructure ( columnsBlockStructure )
{

}

MatrixBlockStructure::MatrixBlockStructure ( const VectorBlockStructure& vectorBlockStructure )
    : M_rowsBlockStructure ( vectorBlockStructure ), M_columnsBlockStructure ( vectorBlockStructure )
{

}

MatrixBlockStructure::MatrixBlockStructure ( const MatrixBlockStructure& blockStructure )
    : M_rowsBlockStructure ( blockStructure.M_rowsBlockStructure ),
      M_columnsBlockStructure ( blockStructure.M_columnsBlockStructure )
{

}

MatrixBlockStructure::~MatrixBlockStructure()
{

}

// ===================================================
// Set Methods
// ===================================================

void
MatrixBlockStructure::setBlockStructure ( const std::vector<UInt>& blockNumRows,
                                          const std::vector<UInt>& blockNumColumns )
{
    M_rowsBlockStructure.setBlockStructure ( blockNumRows );
    M_columnsBlockStructure.setBlockStructure ( blockNumColumns );
}

void
MatrixBlockStructure::setBlockStructure ( const std::vector<UInt>& blocksSize )
{
    M_rowsBlockStructure.setBlockStructure ( blocksSize );
    M_columnsBlockStructure.setBlockStructure ( blocksSize );
}

void
MatrixBlockStructure::setBlockStructure ( const MapVector<MapEpetra>& rowMapVector,
                                          const MapVector<MapEpetra>& columnMapVector )
{
    M_rowsBlockStructure.setBlockStructure ( rowMapVector );
    M_columnsBlockStructure.setBlockStructure ( columnMapVector );
}

void
MatrixBlockStructure::setBlockStructure ( const MapVector<MapEpetra>& mapVector )
{
    M_rowsBlockStructure.setBlockStructure ( mapVector );
    M_columnsBlockStructure.setBlockStructure ( mapVector );
}

void
MatrixBlockStructure::setBlockStructure ( const MatrixBlockStructure& blockStructure )
{
    M_rowsBlockStructure.setBlockStructure ( blockStructure.M_rowsBlockStructure );
    M_columnsBlockStructure.setBlockStructure ( blockStructure.M_columnsBlockStructure );
}

void
MatrixBlockStructure::setBlockStructure ( const VectorBlockStructure& rowsBlockStructure,
                                          const VectorBlockStructure& columnsBlockStructure )
{
    M_rowsBlockStructure.setBlockStructure ( rowsBlockStructure );
    M_columnsBlockStructure.setBlockStructure ( columnsBlockStructure );
}

void
MatrixBlockStructure::setBlockStructure ( const VectorBlockStructure& vectorBlockStructure )
{
    M_rowsBlockStructure.setBlockStructure ( vectorBlockStructure );
    M_columnsBlockStructure.setBlockStructure ( vectorBlockStructure );
}

// ===================================================
// Get Methods
// ===================================================

UInt
MatrixBlockStructure::blockNumRows ( const UInt& rowIndex ) const
{
    return M_rowsBlockStructure.blockSize ( rowIndex );
}

UInt
MatrixBlockStructure::blockNumColumns ( const UInt& columnIndex ) const
{
    return M_columnsBlockStructure.blockSize ( columnIndex );
}

UInt
MatrixBlockStructure::rowBlockFirstIndex ( const UInt& index ) const
{
    return M_rowsBlockStructure.blockFirstIndex ( index );
}

UInt
MatrixBlockStructure::columnBlockFirstIndex ( const UInt& index ) const
{
    return M_columnsBlockStructure.blockFirstIndex ( index );
}

UInt
MatrixBlockStructure::numRowBlocks() const
{
    return M_rowsBlockStructure.numBlocks();
}

UInt
MatrixBlockStructure::numColumnBlocks() const
{
    return M_columnsBlockStructure.numBlocks();
}

UInt
MatrixBlockStructure::numRows() const
{
    return M_rowsBlockStructure.totalSize();
}

UInt
MatrixBlockStructure::numColumns() const
{
    return M_columnsBlockStructure.totalSize();
}

const VectorBlockStructure&
MatrixBlockStructure::rowsBlockStructure() const
{
    return M_rowsBlockStructure;
}

const VectorBlockStructure&
MatrixBlockStructure::columnsBlockStructure() const
{
    return M_columnsBlockStructure;
}

} // Namespace LifeV
