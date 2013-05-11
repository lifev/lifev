//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file

     @brief This file contains the complete specialization of the ETCurrentFE methods updateDetJacobian and updateInverseJacobian.

     @date 12/2012
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
     @author Luca Pasquale <luca.pasquale@mail.polimi.it>
 */

#include <lifev/eta/fem/ETCurrentFE.hpp>

namespace LifeV
{
// ---------------------------------------------------------------
// 1D field
// ---------------------------------------------------------------

// Full specialization for the computation of the determinant
template<>
void
ETCurrentFE<1, 1>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] = M_jacobian[iQuadPt][0][0];
}

// Full specialization for the computation of the determinant
template<>
void
ETCurrentFE<2, 1>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] = M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1]
                             - M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][0][1];
}

// Full specialization for the computation of the determinant
template<>
void
ETCurrentFE<3, 1>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] =
        M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][2]
        + M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][0]
        + M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][1]
        - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][1]
        - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][2]
        - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][0];
}

// Full specialization for the computation of the inverse jacobian
template<>
void
ETCurrentFE<1, 1>::
updateInverseJacobian (const UInt& iQuadPt)
{

    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    M_tInverseJacobian[iQuadPt][0][0] = 1.0 / M_jacobian[iQuadPt][0][0];
}

template<>
void
ETCurrentFE<2, 1>::
updateInverseJacobian (const UInt& iQuadPt)
{

    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    Real det = M_detJacobian[iQuadPt];

    M_tInverseJacobian[iQuadPt][0][0] =  M_jacobian[iQuadPt][1][1] / det;
    M_tInverseJacobian[iQuadPt][1][0] = -M_jacobian[iQuadPt][0][1] / det;
    M_tInverseJacobian[iQuadPt][0][1] = -M_jacobian[iQuadPt][1][0] / det;
    M_tInverseJacobian[iQuadPt][1][1] =  M_jacobian[iQuadPt][0][0] / det;
}

template<>
void
ETCurrentFE<3, 1>::
updateInverseJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    Real det = M_detJacobian[iQuadPt];

    M_tInverseJacobian[iQuadPt][0][0] = ( M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][2]
                                          - M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][1]) / det;

    M_tInverseJacobian[iQuadPt][0][1] = ( M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][0]
                                          - M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][2]) / det;

    M_tInverseJacobian[iQuadPt][0][2] = ( M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][1]
                                          - M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][0]) / det;

    M_tInverseJacobian[iQuadPt][1][0] = ( M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][2][1]
                                          - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][2][2]) / det;

    M_tInverseJacobian[iQuadPt][1][1] = ( M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][2][2]
                                          - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][2][0]) / det;

    M_tInverseJacobian[iQuadPt][1][2] = ( M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][2][0]
                                          - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][2][1]) / det;

    M_tInverseJacobian[iQuadPt][2][0] = ( M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][2]
                                          - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][1]) / det;

    M_tInverseJacobian[iQuadPt][2][1] = ( M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][0]
                                          - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][2]) / det;

    M_tInverseJacobian[iQuadPt][2][2] = ( M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1]
                                          - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][0]) / det;
}

// ---------------------------------------------------------------
// 2D field
// ---------------------------------------------------------------

// Full specialization for the computation of the determinant
template<>
void
ETCurrentFE<1, 2>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] = M_jacobian[iQuadPt][0][0];
}

template<>
void
ETCurrentFE<2, 2>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] = M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1]
                             - M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][0][1];
}

template<>
void
ETCurrentFE<3, 2>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] =
        M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][2]
        + M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][0]
        + M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][1]
        - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][1]
        - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][2]
        - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][0];
}

// Full specialization for the computation of the inverse jacobian
template<>
void
ETCurrentFE<1, 2>::
updateInverseJacobian (const UInt& iQuadPt)
{

    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    M_tInverseJacobian[iQuadPt][0][0] = 1.0 / M_jacobian[iQuadPt][0][0];
}

template<>
void
ETCurrentFE<2, 2>::
updateInverseJacobian (const UInt& iQuadPt)
{

    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    Real det = M_detJacobian[iQuadPt];

    M_tInverseJacobian[iQuadPt][0][0] =  M_jacobian[iQuadPt][1][1] / det;
    M_tInverseJacobian[iQuadPt][1][0] = -M_jacobian[iQuadPt][0][1] / det;
    M_tInverseJacobian[iQuadPt][0][1] = -M_jacobian[iQuadPt][1][0] / det;
    M_tInverseJacobian[iQuadPt][1][1] =  M_jacobian[iQuadPt][0][0] / det;
}

template<>
void
ETCurrentFE<3, 2>::
updateInverseJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    Real det = M_detJacobian[iQuadPt];

    M_tInverseJacobian[iQuadPt][0][0] = ( M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][2]
                                          - M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][1]) / det;

    M_tInverseJacobian[iQuadPt][0][1] = ( M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][0]
                                          - M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][2]) / det;

    M_tInverseJacobian[iQuadPt][0][2] = ( M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][1]
                                          - M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][0]) / det;

    M_tInverseJacobian[iQuadPt][1][0] = ( M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][2][1]
                                          - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][2][2]) / det;

    M_tInverseJacobian[iQuadPt][1][1] = ( M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][2][2]
                                          - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][2][0]) / det;

    M_tInverseJacobian[iQuadPt][1][2] = ( M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][2][0]
                                          - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][2][1]) / det;

    M_tInverseJacobian[iQuadPt][2][0] = ( M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][2]
                                          - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][1]) / det;

    M_tInverseJacobian[iQuadPt][2][1] = ( M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][0]
                                          - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][2]) / det;

    M_tInverseJacobian[iQuadPt][2][2] = ( M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1]
                                          - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][0]) / det;
}

// ---------------------------------------------------------------
// 3D field
// ---------------------------------------------------------------

// Full specialization for the computation of the determinant
template<>
void
ETCurrentFE<1, 3>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] = M_jacobian[iQuadPt][0][0];
}

template<>
void
ETCurrentFE<2, 3>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] = M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1]
                             - M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][0][1];
}

template<>
void
ETCurrentFE<3, 3>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] =
        M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][2]
        + M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][0]
        + M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][1]
        - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][1]
        - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][2]
        - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][0];
}

// Full specialization for the computation of the inverse jacobian
template<>
void
ETCurrentFE<1, 3>::
updateInverseJacobian (const UInt& iQuadPt)
{

    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    M_tInverseJacobian[iQuadPt][0][0] = 1.0 / M_jacobian[iQuadPt][0][0];
}

template<>
void
ETCurrentFE<2, 3>::
updateInverseJacobian (const UInt& iQuadPt)
{

    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    Real det = M_detJacobian[iQuadPt];

    M_tInverseJacobian[iQuadPt][0][0] =  M_jacobian[iQuadPt][1][1] / det;
    M_tInverseJacobian[iQuadPt][1][0] = -M_jacobian[iQuadPt][0][1] / det;
    M_tInverseJacobian[iQuadPt][0][1] = -M_jacobian[iQuadPt][1][0] / det;
    M_tInverseJacobian[iQuadPt][1][1] =  M_jacobian[iQuadPt][0][0] / det;
}

template<>
void
ETCurrentFE<3, 3>::
updateInverseJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    Real det = M_detJacobian[iQuadPt];

    M_tInverseJacobian[iQuadPt][0][0] = ( M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][2]
                                          - M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][1]) / det;

    M_tInverseJacobian[iQuadPt][0][1] = ( M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][0]
                                          - M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][2]) / det;

    M_tInverseJacobian[iQuadPt][0][2] = ( M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][1]
                                          - M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][0]) / det;

    M_tInverseJacobian[iQuadPt][1][0] = ( M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][2][1]
                                          - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][2][2]) / det;

    M_tInverseJacobian[iQuadPt][1][1] = ( M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][2][2]
                                          - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][2][0]) / det;

    M_tInverseJacobian[iQuadPt][1][2] = ( M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][2][0]
                                          - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][2][1]) / det;

    M_tInverseJacobian[iQuadPt][2][0] = ( M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][2]
                                          - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][1]) / det;

    M_tInverseJacobian[iQuadPt][2][1] = ( M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][0]
                                          - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][2]) / det;

    M_tInverseJacobian[iQuadPt][2][2] = ( M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1]
                                          - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][0]) / det;
}


} // Namespace LifeV
