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

#include <life/lifemesh/geo0D.hpp>
#include <lifeconfig.h>

namespace LifeV
{


// ==========================================
// Constructor & Destructor
// ==========================================
Geo0D::Geo0D() :
		MeshEntityWithBoundary( 0 ),
        M_coordinates()
{
    M_coordinates.assign( 0 );
}

Geo0D::Geo0D( ID identity, bool boundary )
        :
        MeshEntityWithBoundary( identity, boundary ),
        M_coordinates()
{
    M_coordinates.assign( 0 );
}

Geo0D::Geo0D( ID identity, Real x, Real y, Real z, bool boundary )
        :
        MeshEntityWithBoundary( identity, boundary ),
        M_coordinates()
{
	M_coordinates[ 0 ] = x;
	M_coordinates[ 1 ] = y;
	M_coordinates[ 2 ] = z;
}

Geo0D::Geo0D( Geo0D const & Element )
        :
        MeshEntityWithBoundary( Element.id(), Element.boundary() ),
        M_coordinates( Element.M_coordinates )
{

}

// ==========================================
// Operators
// ==========================================
Geo0D &
Geo0D::operator=( Geo0D const & Element )
{
    if (  this == &Element )
        return *this;

    this->setId(Element.id());
    setBoundary(Element.boundary());
    M_coordinates = Element.M_coordinates;
    return *this;
}

// ==========================================
// Methods
// ==========================================
std::ostream &
Geo0D::showMe( bool verbose, std::ostream & out ) const
{
    out.setf( std::ios::scientific, std::ios::floatfield );
    out << "----- Geo0D object -----" << std::endl;
    if ( verbose )
    {
        unsigned i;
        out << " Coordinates:" << std::endl;
        Real const * coordinateVector = coor();
        for ( i = 0; i < nDimensions-1; i++ )
        {
            out << coordinateVector[ i ] << ",  ";
        }
        out << coordinateVector[ i ] << std::endl << std::endl;
    }
    out << " ID       = " << id()      << std::endl;
    out << " local ID = " << localId() << std::endl;

    out << "----- END OF Geo0D data ---" << std::endl << std::endl;
    return out;
}

}
