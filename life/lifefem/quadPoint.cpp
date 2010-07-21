//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER


#include <quadPoint.hpp>

namespace LifeV {

QuadPoint::QuadPoint()
    : M_weight( 0 ), 
      M_coor( 3 )
{
    M_coor[ 0 ] = 0;
    M_coor[ 1 ] = 0;
    M_coor[ 2 ] = 0;
}

QuadPoint::QuadPoint( Real x, Real y, Real z, Real weight )
    : M_weight( weight ),
      M_coor( 3 )
{
    M_coor[ 0 ] = x;
    M_coor[ 1 ] = y;
    M_coor[ 2 ] = z;
}

QuadPoint::QuadPoint( Real x, Real y, Real weight )
    : M_weight( weight ),
      M_coor( 3 )
{
    M_coor[ 0 ] = x;
    M_coor[ 1 ] = y;
    M_coor[ 2 ] = 0.;
}

QuadPoint::QuadPoint( Real x, Real weight )
    : M_weight( weight ),
      M_coor( 3 )
{
    M_coor[ 0 ] = x;
    M_coor[ 1 ] = 0.;
    M_coor[ 2 ] = 0.;
}

QuadPoint::QuadPoint(const GeoVector& coor, const Real& weight)
    : M_weight(weight),
      M_coor(coor)
{}

QuadPoint::QuadPoint(const GeoVector& coor, const Real& weight, const UInt& spaceDim)
    : M_weight(weight),
      M_coor(spaceDim)
{
    for (UInt i(0); (i<spaceDim) && (i<coor.size()); ++i)
    {
        M_coor[i] = coor[i];
    }
    
    // Add zeros if necessary
    for (UInt i(coor.size()); i<spaceDim ; ++i)
    {
        M_coor[i] = 0.0;
    }
}

QuadPoint::QuadPoint(const QuadPoint& qp)
    : M_weight(qp.M_weight),
      M_coor(qp.M_coor.size())
{
    M_coor = qp.M_coor;
}

QuadPoint::QuadPoint(const QuadPoint& qp, const UInt spaceDim)
    : M_weight(qp.M_weight),
      M_coor(spaceDim)
{
    for (UInt i(0); (i<spaceDim) && (i<qp.M_coor.size()); ++i)
    {
        M_coor[i] = qp.M_coor[i];
    }
    
    // Add zeros if necessary
    for (UInt i(qp.M_coor.size()); i<spaceDim ; ++i)
    {
        M_coor[i] = 0.0;
    }
}

} // Namespace LifeV
