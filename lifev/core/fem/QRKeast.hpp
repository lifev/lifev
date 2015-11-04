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
    @brief File for the definition of the quadrature rule defined in
    the article [Keast85]

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 12-2011

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#ifndef QR_KEAST_HPP
#define QR_KEAST_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/core/mesh/ElementShapes.hpp>

#include <lifev/core/fem/QuadraturePoint.hpp>

#include <iostream>
#include <fstream>

namespace LifeV
{


//! QRKeast - A set of quadrature for tetrahedra
/*!
  We concentrate on the lowest degree, so negative weights
  are possible!
*/

template< UInt degree >
class QRKeast
{
public:

private:
    QRKeast();
    ~QRKeast();
    QRKeast ( const QRKeast&);
};


/* Specializations */
template <>
class QRKeast<1>
{
public:

    QRKeast() {};
    ~QRKeast() {};

    static const Real& weight (const UInt& iPt)
    {
        return M_weights[iPt];
    }

    static const Real& quadPointCoor ( const UInt iPt, const UInt iCoor)
    {
        return M_points[iPt][iCoor];
    }

    static VectorSmall<3> quadPointCoor ( const UInt iPt)
    {
        return VectorSmall<3> (M_points[iPt][0], M_points[iPt][1], M_points[iPt][2]);
    }

    static QuadraturePoint quadPoint (const UInt iPt)
    {
        return QuadraturePoint (M_points[iPt][0], M_points[iPt][1], M_points[iPt][2], M_weights[iPt]);
    }

    static UInt nbQuadPt()
    {
        return 1;
    }

    static ReferenceShapes shape()
    {
        return TETRA;
    }

    static UInt dimension()
    {
        return 3;
    }

private:

    static const Real M_points[1][3];
    static const Real M_weights[1];
};

template <>
class QRKeast<4>
{
public:

    QRKeast() {};
    ~QRKeast() {};

    static const Real& weight (const UInt& iPt)
    {
        return M_weights[iPt];
    }

    static const Real& quadPointCoor ( const UInt iPt, const UInt iCoor)
    {
        return M_points[iPt][iCoor];
    }

    static VectorSmall<3> quadPointCoor ( const UInt iPt)
    {
        return VectorSmall<3> (M_points[iPt][0], M_points[iPt][1], M_points[iPt][2]);
    }

    static QuadraturePoint quadPoint (const UInt iPt)
    {
        return QuadraturePoint (M_points[iPt][0], M_points[iPt][1], M_points[iPt][2], M_weights[iPt]);
    }

    static UInt nbQuadPt()
    {
        return 11;
    }

    static ReferenceShapes shape()
    {
        return TETRA;
    }

    static UInt dimension()
    {
        return 3;
    }

private:

    static const Real M_points[11][3];
    static const Real M_weights[11];
};

template <>
class QRKeast<6>
{
public:

    QRKeast() {};
    ~QRKeast() {};

    static const Real& weight (const UInt& iPt)
    {
        return M_weights[iPt];
    }

    static const Real& quadPointCoor ( const UInt iPt, const UInt iCoor)
    {
        return M_points[iPt][iCoor];
    }

    static VectorSmall<3> quadPointCoor ( const UInt iPt)
    {
        return VectorSmall<3> (M_points[iPt][0], M_points[iPt][1], M_points[iPt][2]);
    }

    static QuadraturePoint quadPoint (const UInt iPt)
    {
        return QuadraturePoint (M_points[iPt][0], M_points[iPt][1], M_points[iPt][2], M_weights[iPt]);
    }

    static UInt nbQuadPt()
    {
        return 24;
    }

    static ReferenceShapes shape()
    {
        return TETRA;
    }

    static UInt dimension()
    {
        return 3;
    }

private:

    static const Real M_points[24][3];
    static const Real M_weights[24];
};


template <>
class QRKeast<7>
{
public:

    QRKeast() {};
    ~QRKeast() {};

    static const Real& weight (const UInt& iPt)
    {
        return M_weights[iPt];
    }

    static const Real& quadPointCoor ( const UInt iPt, const UInt iCoor)
    {
        return M_points[iPt][iCoor];
    }

    static VectorSmall<3> quadPointCoor ( const UInt iPt)
    {
        return VectorSmall<3> (M_points[iPt][0], M_points[iPt][1], M_points[iPt][2]);
    }

    static QuadraturePoint quadPoint (const UInt iPt)
    {
        return QuadraturePoint (M_points[iPt][0], M_points[iPt][1], M_points[iPt][2], M_weights[iPt]);
    }

    static UInt nbQuadPt()
    {
        return 31;
    }

    static ReferenceShapes shape()
    {
        return TETRA;
    }

    static UInt dimension()
    {
        return 3;
    }

private:

    static const Real M_points[31][3];
    static const Real M_weights[31];
};


}
#endif
