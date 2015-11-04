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
    @brief Zero dimensional entity

    @author Luca Formaggia <luca.formaggia@polimi.it>
    @contributor Marta D'Elia <mdelia2@mathcs.emory.edu>
    @maintainer Marta D'Elia <mdelia2@mathcs.emory.edu>

    @date 00-00-0000

*/

#include <lifev/core/mesh/MeshVertex.hpp>
#include <lifev/core/LifeV.hpp>

namespace LifeV
{


// ==========================================
// Constructor & Destructor
// ==========================================
MeshVertex::MeshVertex() : M_coordinates()
{
}

MeshVertex::MeshVertex ( ID identity, bool boundary )
    : MeshEntity ( identity, boundary ), M_coordinates()
{
}

MeshVertex::MeshVertex ( ID identity, Real x, Real y, Real z, bool boundary )
    : MeshEntity ( identity, boundary ), M_coordinates ( x, y, z )
{
}

// ==========================================
// Operators
// ==========================================


// ==========================================
// Methods
// ==========================================
std::ostream&
MeshVertex::showMe ( bool verbose, std::ostream& out ) const
{
    out.setf ( std::ios::scientific, std::ios::floatfield );
    out << "----- MeshVertex object -----" << std::endl;
    if ( verbose )
    {
        UInt i;
        out << " Coordinates:" << std::endl;
        Real const* coordinateVector = coordinatesArray();
        for ( i = 0; i < nDimensions - 1; i++ )
        {
            out << coordinateVector[ i ] << ",  ";
        }
        out << coordinateVector[ i ] << std::endl << std::endl;
    }
    MeshEntity::showMe ( out );

    out << "----- END OF MeshVertex data ---" << std::endl << std::endl;
    return out;
}

}
