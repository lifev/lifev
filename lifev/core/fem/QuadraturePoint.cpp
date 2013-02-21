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
    @brief Implementation of the quadPoint class, usefull for the quadrature rules.

    @author Jean-Frederic Gerbeau
            Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 18-05-2010

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#include <lifev/core/fem/QuadraturePoint.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

QuadraturePoint::QuadraturePoint()
    : M_weight ( 0 ),
      M_coor ( 3 )
{
    M_coor[ 0 ] = 0;
    M_coor[ 1 ] = 0;
    M_coor[ 2 ] = 0;
}

QuadraturePoint::QuadraturePoint ( Real x, Real y, Real z, Real weight )
    : M_weight ( weight ),
      M_coor ( 3 )
{
    M_coor[ 0 ] = x;
    M_coor[ 1 ] = y;
    M_coor[ 2 ] = z;
}

QuadraturePoint::QuadraturePoint ( Real x, Real y, Real weight )
    : M_weight ( weight ),
      M_coor ( 3 )
{
    M_coor[ 0 ] = x;
    M_coor[ 1 ] = y;
    M_coor[ 2 ] = 0.;
}

QuadraturePoint::QuadraturePoint ( Real x, Real weight )
    : M_weight ( weight ),
      M_coor ( 3 )
{
    M_coor[ 0 ] = x;
    M_coor[ 1 ] = 0.;
    M_coor[ 2 ] = 0.;
}

QuadraturePoint::QuadraturePoint (const GeoVector& coor, const Real& weight)
    : M_weight (weight),
      M_coor (coor)
{}

QuadraturePoint::QuadraturePoint (const GeoVector& coor, const Real& weight, const UInt& spaceDim)
    : M_weight (weight),
      M_coor (spaceDim)
{
    for (UInt i (0); (i < spaceDim) && (i < coor.size() ); ++i)
    {
        M_coor[i] = coor[i];
    }

    // Add zeros if necessary
    for (UInt i (coor.size() ); i < spaceDim ; ++i)
    {
        M_coor[i] = 0.0;
    }
}

QuadraturePoint::QuadraturePoint (const QuadraturePoint& qp)
    : M_weight (qp.M_weight),
      M_coor (qp.M_coor.size() )
{
    M_coor = qp.M_coor;
}

QuadraturePoint::QuadraturePoint (const QuadraturePoint& qp, const UInt spaceDim)
    : M_weight (qp.M_weight),
      M_coor (spaceDim)
{
    for (UInt i (0); (i < spaceDim) && (i < qp.M_coor.size() ); ++i)
    {
        M_coor[i] = qp.M_coor[i];
    }

    // Add zeros if necessary
    for (UInt i (qp.M_coor.size() ); i < spaceDim ; ++i)
    {
        M_coor[i] = 0.0;
    }
}

} // Namespace LifeV
