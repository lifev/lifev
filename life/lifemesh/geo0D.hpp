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
/*! file geo0D.h */
#ifndef _GEO0D_HH_
#define _GEO0D_HH_

#include <boost/array.hpp>

#include "meshEntity.hpp"
#include "basisElSh.hpp"

namespace LifeV
{
//! \defgroup GeoXD Basis Geometrical Entities Geo0D and GeoND.
/*!

They are intermediate classes used to build the actual Geometry classes

\warning Geo1D/2D/3D are template classes because some of the info is not
known a priori and I want all vector dimensions determined at compile time
to enhance memory access time. IT IS UP TO THE USER to provide a coherent
GeoShape!.  */

/*@{*/

//! Zero dimensional entity. It stores boundary information
class Geo0D
    :
        public MeshEntityWithBoundary
{
public:

    Geo0D();
    //! constructor where I give the id and declare if Geo0D object is on a
    //!boundary
    explicit Geo0D( ID id, bool boundary = false );
    //! constructor where I give the id, the point coordinate and I declare
    //! if the Geo0D object is on a boundary
    explicit Geo0D( ID id, Real x, Real y, Real z, bool boundary = false );

    Geo0D( Geo0D const & G );
    Geo0D & operator=( Geo0D const & G );

    typedef GeoPoint GeoShape;


    //! returns a pointer to a Real[3] containing the coordinates
    Real * coor()
        {
            return _coor.c_array();
        };
    Real const * coor() const
        {
            return _coor.data();
        };


    //! Used to provide coords to object created using
    //! a constructor with no coordinates given, or to
    //! modify existing coordinates.
    Real & x()
        {
            return _coor[ 0 ];
        }
    Real & y()
        {
            return _coor[ 1 ];
        }
    Real & z()
        {
#if defined(THREEDIM)
            return _coor[ 2 ];
#else

            ERROR_MSG( "z coordinate may be modified only in a 3D problem" );
#endif

        }
    Real x() const
        {
            return _coor[ 0 ];
        }
    Real y() const
        {
            return _coor[ 1 ];
        };
    Real z() const
        {
#if defined(THREEDIM)
            return _coor[ 2 ];
#else

            return 0;
#endif

        }

    //!Another way to access coordinate data
    Real coordinate ( ID const i ) const
        {
            ASSERT_BD( i > 0 && i <= NDIM ) ;
            return _coor[ i -1 ]; // indexing from 1
        }
    //!Another way to modify coordinate data
    Real & coordinate ( ID const i )
        {
            ASSERT_BD( i > 0 && i <= NDIM ) ;
            return _coor[ i -1 ];
        }
    //! Useful for debugging
    std::ostream & showMe( bool verbose = false, std::ostream & c = std::cout ) const;

private:
    boost::array<Real,nDimensions> _coor;
};
}
#endif
