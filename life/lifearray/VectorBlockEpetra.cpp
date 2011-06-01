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
VectorBlockEpetra( const map_type& map, const mapType_type& mapType = Unique)
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
}


// ===================================================
// Operators
// ===================================================

// ===================================================
// Methods
// ===================================================

// ===================================================
// Set Methods
// ===================================================

// ===================================================
// Get Methods
// ===================================================

// ===================================================
// Private Methods
// ===================================================

} // Namespace LifeV
