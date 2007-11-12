/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <life/lifemesh/geo0D.hpp>

namespace LifeV
{
/*--------------------------------------------------------------
                                 Geo0D
  ---------------------------------------------------------------*/
Geo0D::Geo0D()
    :
    MeshEntityWithBoundary( 0 ),
    _coor()
{
    _coor.assign( 0 );
}

Geo0D::Geo0D( ID id, bool boundary )
    :
    MeshEntityWithBoundary( id, boundary ),
    _coor()
{
    _coor.assign( 0 );
}

Geo0D::Geo0D( ID id, Real x, Real y, Real z, bool boundary )
    :
    MeshEntityWithBoundary( id, boundary ),
    _coor()
{
    _coor[ 0 ] = x;
    _coor[ 1 ] = y;
    _coor[ 2 ] = z;
}

Geo0D::Geo0D( Geo0D const & G )
    :
    MeshEntityWithBoundary( G.id(), G._boundary ),
    _coor( G._coor )
{
}

Geo0D &
Geo0D::operator=( Geo0D const & G )
{
    if (  this == &G )
        return *this;
    this->setId(G.id());
    _boundary = G._boundary;
    _coor = G._coor;
    return *this;
}


std::ostream &
Geo0D::showMe( bool verbose, std::ostream & out ) const
{
    out.setf( std::ios::scientific, std::ios::floatfield );
    out << " Geo0D object " << std::endl;
    if ( verbose )
    {
        unsigned i;
        out << " Coordinates:" << std::endl;
        Real const * c = coor();
        for ( i = 0; i < nDimensions-1; i++ )
        {
            out << c[ i ] << ",  ";
        }
        out << c[i] << std::endl << std::endl;
    }
    out << " ID       = " << id()      << std::endl;
    out << " local ID = " << localId() << std::endl;

    out << "----- END OF Geo0D data ---" << std::endl << std::endl;
    return out;
}

}
