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
    @brief  This file contains the the VectorEpetraStructuredView implementation.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 07 Jun 2011
 */

#include <lifev/core/array/VectorEpetraStructuredView.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

VectorEpetraStructuredView::
VectorEpetraStructuredView()
    : M_blockSize(),
      M_firstIndex(),
      M_lastValidIndex(),
      M_vector()
{}

VectorEpetraStructuredView::
VectorEpetraStructuredView ( const VectorEpetraStructuredView& otherView)
    : M_blockSize (otherView.M_blockSize),
      M_firstIndex (otherView.M_firstIndex),
      M_lastValidIndex (otherView.M_lastValidIndex),
      M_vector (otherView.M_vector)
{}

VectorEpetraStructuredView::
~VectorEpetraStructuredView()
{}

// ===================================================
// Methods
// ===================================================

void
VectorEpetraStructuredView::
showMe ( std::ostream& output ) const
{
    output << "VectorBlockViewEpetra informations:" << std::endl
           << "Size = " << M_blockSize << std::endl
           << "First index = " << M_firstIndex << std::endl
           << "Last (valid) index = " << M_lastValidIndex << std::endl;
}

Int
VectorEpetraStructuredView::
sumIntoGlobalValues ( const Int GID, const Real value ) const
{
    ASSERT (GID < static_cast<Int> (M_blockSize), " Error in assembling the block vector: global id to large for the block")

    // Compute the global ID in the monolithic vector:
    // size of the block + location in the block
    const Int TotalGID (GID + M_firstIndex);
    return M_vector->sumIntoGlobalValues ( TotalGID, value );
}

// ===================================================
// Set Methods
// ===================================================

void
VectorEpetraStructuredView::
setup ( const UInt& firstIndex, const UInt& blockSize, vector_Type* vector )
{
    M_blockSize = blockSize;
    M_firstIndex = firstIndex;
    M_lastValidIndex = firstIndex + blockSize - 1;
    M_vector = vector;
}

} // Namespace LifeV
