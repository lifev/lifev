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
    @brief File containing the procedures for the local assembly of the differential operators for the structural problem

    @author Gianmarco Mengaldo <gianmarco.mengaldo@gmail.com>
    @author Paolo Tricerri <gianmarco.mengaldo@gmail.com>
    @mantainer Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef ELEMOPERSTRUCTURE_CPP
#define ELEMOPERSTRUCTURE_CPP 1

#include <lifev/structure/fem/AssemblyElementalStructure.hpp>
#include <boost/multi_array.hpp>

namespace LifeV
{

namespace AssemblyElementalStructure
{

void computeGradientLocalDisplacement (boost::multi_array<Real, 3>& gradientLocalDisplacement,
                                       const VectorElemental& uk_loc, const CurrentFE& fe )
{
    // \grad u^k at each quadrature poInt
    Real s;

    // loop on quadrature poInts
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < nDimensions; icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < nDimensions; jcoor++ )
            {
                s = 0.0;
                for (UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    //  \grad u^k at a quadrature poInt
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
                }
                gradientLocalDisplacement[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
}

void computeLocalDeformationGradient (const VectorElemental& uk_loc, std::vector<Epetra_SerialDenseMatrix>& tensorF, const CurrentFE& fe )
{
    // \grad u^k at each quadrature poInt
    Real s;

    for ( Int k=0; k < static_cast<Int> (fe.nbQuadPt()); k++)
    {
        // loop on space coordinates
        for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); icoor++ )
        {
            // loop  on space coordinates
            for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); jcoor++ )
            {
                s = 0.0;
                for (Int i = 0; i < static_cast<Int> (fe.nbFEDof()); i++ )
                {
                    //  \grad u^k at a quadrature poInt
                    s += fe.phiDer( i, jcoor, k ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
                }
                tensorF[k]( icoor, jcoor ) = s;

                if ( icoor == jcoor )
                    tensorF[k] ( icoor, jcoor ) += 1.0;
            }
        }
    }
}

void computeLocalDeformationGradientWithoutIdentity (const VectorElemental& uk_loc, std::vector<Epetra_SerialDenseMatrix>& tensorF, const CurrentFE& fe )
{
    // \grad u^k at each quadrature poInt
    Real s;

    for ( Int k=0; k < static_cast<Int> (fe.nbQuadPt()); k++)
    {
        // loop on space coordinates
        for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); icoor++ )
        {
            // loop  on space coordinates
            for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); jcoor++ )
            {
                s = 0.0;
                for (Int i = 0; i < static_cast<Int> (fe.nbFEDof()); i++ )
                {
                    //  \grad u^k at a quadrature poInt
                    s += fe.phiDer( i, jcoor, k ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
                }
                tensorF[k]( icoor, jcoor ) = s;
            }
        }
    }
}

// The methods for linear elastic model (stiff_strain and stiff_div) are implemented in AssemblyElemental.cpp

//! Methods for St. Venant Kirchhoff model
//! Methods for the stiffness matrix

//! \f$ coef \cdot ( trace { [\nabla u^k]^T \nabla u }, \nabla\cdot  v  ) \f$
void stiff_derdiv( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe )
{
    Real s;

    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
                        {
                            s += fe.phiDer ( i, icoor, ig ) * gradientLocalDisplacement[ jcoor ][ k ][ ig ]
                                 * fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                        }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}



//! \f$ coef \cdot ( [\nabla u^k]^T \nabla u : \nabla v  )\f$
void stiff_dergradbis ( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
                        MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;

    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor <  fe.nbCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            s += fe.phiDer ( i, k, ig ) * gradientLocalDisplacement[ jcoor][ icoor ][ ig ] *
                                 fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}



//! \f$  coef * ( (\div u_k) \grad u : \grad v  )
void stiff_divgrad ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;
    // local vector for \div u^k at each quadrature poInt
    //Real duk[ fe.nbQuadPt() ];
    std::vector<Real > duk (fe.nbQuadPt(), 0.0);

    // loop on quadrature poInts
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        s = 0;
        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {
            for ( UInt i = 0; i < fe.nbFEDof(); i++ )
            {
                // construction of \div u^k at a quadrature poInt
                s += fe.phiDer ( i, icoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
        }
        duk[ ig ] = s;
    }

    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0.0;
            for ( UInt k = 0; k < fe.nbCoor(); ++k )
            {
                for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += duk[ ig ] * fe.phiDer ( i, k, ig ) *  fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }

}

//! \f$ coef * ( \grad u_k : \grad u_k) * ( \grad u : \grad v  )
void stiff_gradgrad ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{
    Real s,s1;
    //    (\grad u_k : \grad u_k) at each quadrature poInt
    //Real gguk[ fe.nbQuadPt() ];
    std::vector<Real > gguk (fe.nbQuadPt(), 0.0);

    // loop on quadrature poInts
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        s = 0;
        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {
            for ( UInt l = 0; l < fe.nbCoor(); l++ )
            {
                s1 = 0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s1 += fe.phiDer ( i, l, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
                }
                s += s1 * s1;
            }
        }
        gguk[ ig ] = s;
    }

    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0.0;
            for ( UInt k = 0; k < fe.nbCoor(); ++k )
            {
                for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += gguk[ ig ] * fe.phiDer ( i, k, ig ) *  fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }
}

void stiff_dergrad_gradbis ( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
                             MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;

    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                    {
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            s += gradientLocalDisplacement[ icoor ][ jcoor ][ ig ] * fe.phiDer ( i, k, ig ) *
                                 fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}



// coef * ( \grad u^k [\grad u]^T : \grad v )
void stiff_dergrad_gradbis_Tr ( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
                                MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;
    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                    {
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            s += gradientLocalDisplacement[ icoor ][ k ][ ig ]  * fe.phiDer ( j, k, ig ) *
                                 fe.phiDer ( i, jcoor, ig ) * fe.weightDet ( ig );
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}



// coef * ( \grad u^k [\grad u^k]^T \grad u : \grad v )
void stiff_gradgradTr_gradbis ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;

    //! \grad u^k  [\grad u^k]^T  at each quadrature poInt
    //Real guk_gukT[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];
    boost::multi_array<Real, 3> guk_gukT (boost::extents[fe.nbCoor()][fe.nbCoor()][fe.nbQuadPt()]);

    // loop on quadrature poInts
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt n = 0; n < fe.nbCoor(); n++ )
                {
                    for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    {
                        for ( UInt j = 0; j < fe.nbFEDof(); j++ )
                        {
                            //! \grad u^k  [\grad u^k]^T  at each quadrature poInt
                            s  += fe.phiDer ( i, n, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ] *
                                  fe.phiDer ( j, n, ig ) * uk_loc.vec() [ j + jcoor * fe.nbFEDof() ];
                        }
                    }
                }
                guk_gukT[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            s += fe.phiDer ( i, k, ig ) * guk_gukT[ icoor ][ jcoor ][ ig ] *
                                 fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}
// End of methods for the stiffness matrix (St. Venant-Kirchhoff material)

//! Methods for the jacobian (St. Venant-Kirchhoff material)

//! \f$ coef \cdot ( [\nabla u]^T \nabla u^k + [\nabla u^k]^T \nabla u : \nabla v  )\f$
void stiff_dergrad ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement,
                     MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                        {
                            s += fe.phiDer ( i, k, ig ) *
                                 ( gradientLocalDisplacement[ jcoor ][ k ][ ig ] * fe.phiDer ( j, icoor, ig )
                                   + gradientLocalDisplacement[ jcoor ][ icoor ][ ig ] * fe.phiDer ( j, k, ig ) ) * fe.weightDet ( ig );
                        }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}



// coef * ( (\div u) \grad u_k : \grad v  )
void stiff_divgrad_2 ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement,
                       MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;

    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            s += fe.phiDer ( j, jcoor, ig ) *
                                 gradientLocalDisplacement[ icoor ][ k ][ ig ] *
                                 fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}

// coef * ( \grad u_k : \grad u) *( \grad u_k : \grad v  )
void stiff_gradgrad_2( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                    {
                        for ( UInt l = 0; l < fe.nbCoor(); ++l )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                                s += gradientLocalDisplacement[ jcoor ][ l ][ ig ] * fe.phiDer ( j, l, ig ) *
                                     gradientLocalDisplacement[ icoor ][ k ][ ig ] * fe.phiDer (i, k, ig ) *
                                     fe.weightDet ( ig );
                        }
                    }
                    mat ( i, j ) += coef  * s;
                }
            }
        }
    }
}

// coef * ( \grad \delta u [\grad u^k]^T : \grad v )
void stiff_dergrad_gradbis_2 ( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
                               MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;

    //
    // blocks (icoor,jcoor) of elmat
    //
    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0.0;
            for ( UInt l = 0; l < fe.nbCoor(); ++l )
            {
                for ( UInt k = 0; k < fe.nbCoor(); ++k )
                {
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                        s += gradientLocalDisplacement[ l ][ k ][ ig ] * fe.phiDer ( i, k, ig ) *
                             fe.phiDer ( j, l, ig ) * fe.weightDet ( ig );
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }
}

void stiff_dergrad_gradbis_Tr_2 ( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
                                  MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;
    //
    // blocks (icoor,jcoor) of elmat
    //
    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0.0;
            for ( UInt l = 0; l < fe.nbCoor(); ++l )
            {
                for ( UInt k = 0; k < fe.nbCoor(); ++k )
                {
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    {
                        s += gradientLocalDisplacement[ k ][ l ][ ig ] * fe.phiDer ( i, k, ig ) *
                             fe.phiDer ( j, l, ig ) * fe.weightDet ( ig );
                    }
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }
}



// coef * (  \grad u^k [\grad u]^T \grad u^k : \grad v  )
void stiff_gradgradTr_gradbis_2 ( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
                                  MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;

    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor ); // it extracts the (icoor, jcoor) block

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt l = 0; l < fe.nbCoor(); ++l )
                    {
                        for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += gradientLocalDisplacement[ icoor ][ l ][ ig ] *
                                     gradientLocalDisplacement[ jcoor ][ k ][ ig ] * fe.phiDer ( i, k, ig ) *
                                     fe.phiDer ( j, l, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}
//  coef * (  \grad u [\grad u^k]^T \grad u^k : \grad v  )
void stiff_gradgradTr_gradbis_3 ( Real coef, const VectorElemental& uk_loc,
                                  MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;
    // \grad u^k  [\grad u^k]^T  at each quadrature poInt
    //Real guk_gukT[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];
    boost::multi_array<Real, 3> guk_gukT (boost::extents[fe.nbCoor()][fe.nbCoor()][fe.nbQuadPt()]);

    // loop on quadrature poInts                                                // (\grad u^k  [\grad u^k]^T )^T
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt n = 0; n < fe.nbCoor(); n++ )
                {
                    for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    {
                        for ( UInt j = 0; j < fe.nbFEDof(); j++ )
                        {
                            // \grad u^k  [\grad u^k]^T  at each quadrature poInt
                            s  += fe.phiDer ( i, n, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ] *
                                  fe.phiDer ( j, n, ig ) * uk_loc.vec() [ j + jcoor * fe.nbFEDof() ] ;
                        }
                    }
                }
                guk_gukT[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }

    //
    // blocks (icoor,jcoor) of elmat
    //
    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0.0;
            for ( UInt l = 0; l < fe.nbCoor(); ++l )
            {
                for ( UInt k = 0; k < fe.nbCoor(); ++k )
                {
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    {
                        s += guk_gukT[ k ][ l ][ ig ] * fe.phiDer ( i, k, ig ) *  fe.phiDer ( j, l, ig ) * fe.weightDet ( ig );
                    }
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }
}
// End of St. Venant Kirchhoff model

//! ***********************************************************************************************
//! METHODS SHARED BETWEEN NEO-HOOKEAN AND EXPONENTIAL MODELS
//! ***********************************************************************************************
//! The volumetric part is the same for Neo-Hookean and Expoential models is the same

//! STIFFFNESS VECTOR -----------------------------------------------------------------------------
//! Volumetric part--------------------------------------------------------------------------------
//! Source term source_Pvol: Int { coef /2* (J^2 - J + log(J) ) * 1/J * (CofF : \nabla v) }
void source_Pvol ( Real      coef,
                   const boost::multi_array<Real, 3 >& CofFk,
                   const std::vector<Real>&  Jk,
                   VectorElemental&  elvec,
                   const CurrentFE&  fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        // block (icoor) of elvec
        VectorElemental::vector_view vec =  elvec.block ( icoor );
        for ( UInt i = 0; i < fe.nbFEDof(); ++i )
        {
            s = 0.0;
            for ( UInt k = 0; k < nDimensions; ++k )
            {
                for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += ( Jk[ig] * Jk[ig] - Jk[ig] + log ( Jk[ig] ) ) * ( 1 / Jk[ig] ) *
                         CofFk[icoor][ k][ ig] * fe.phiDer (i, k, ig) * fe.weightDet (ig);
                }
            }
            vec (i) += coef * s;
        }
    }
}

//! JACOBIAN MATRIX -------------------------------------------------------------------------------
//! Volumetric part--------------------------------------------------------------------------------

//! 1. Jacobian matrix: Int { 1/2 * coef * ( 2 - 1/J + 1/J^2 ) * ( CofF : \nabla \delta ) (CofF : \nabla v) }
void stiff_Jac_Pvol_1term ( Real          coef,
                            const boost::multi_array<Real, 3 >& CofFk,
                            const std::vector<Real>&      Jk,
                            MatrixElemental& elmat,
                            const CurrentFE& fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt l = 0; l < nDimensions; ++l )
                    {
                        for ( UInt k = 0; k < nDimensions; ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( 2.0 - ( 1 / Jk[ig] ) + ( 1 / ( Jk[ig] * Jk[ig] ) ) ) *
                                     CofFk[ jcoor ][ l ][ ig ] * fe.phiDer ( j, l, ig ) *
                                     CofFk[ icoor ][ k ][ ig ] * fe.phiDer ( i, k, ig ) *
                                     fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += s * coef;
                }
            }
        }
    }
}


//! 2. Stiffness matrix: Int { 1/2 * coef * ( 1/J - 1 - log(J)/J^2 ) * ( CofF [\nabla \delta]^t CofF ) : \nabla v }
void stiff_Jac_Pvol_2term ( Real           coef,
                            const boost::multi_array<Real, 3 >&  CofFk,
                            const std::vector<Real>&       Jk,
                            MatrixElemental&  elmat,
                            const CurrentFE&  fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt l = 0; l < nDimensions; ++l )
                    {
                        for ( UInt k = 0; k < nDimensions; ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( ( 1 / Jk[ig] ) - 1. - ( 1 / ( Jk[ig] * Jk[ig] ) ) * log ( Jk[ig] ) ) *
                                     CofFk[ icoor ][ l ][ ig ] * fe.phiDer ( j, l, ig ) *
                                     CofFk[ jcoor ][ k ][ ig ] * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += s * coef;
                }
            }
        }
    }
}


//! ***********************************************************************************************
//! METHODS FOR NEO-HOOKEAN MODEL
//! ***********************************************************************************************
//! Stiffness vector isochoric part ---------------------------------------------------------------

//! Source term source_P1iso_NH: Int { coef * (  J^(-2/3) * (F : \nabla v) - 1/3 * (Ic_iso / J) (CofF : \nabla v) ) }
void source_P1iso_NH ( Real      coef,
                       const boost::multi_array<Real, 3 >& CofFk,
                       const boost::multi_array<Real, 3 >& Fk,
                       const std::vector<Real>&   Jk,
                       const std::vector<Real>&   Ic_isok ,
                       VectorElemental& elvec,
                       const CurrentFE& fe )
{
    Real s1, s2;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        VectorElemental::vector_view vec =  elvec.block ( icoor );
        for ( UInt i = 0; i < fe.nbFEDof(); ++i )
        {
            s1 = 0.0;
            s2 = 0.0;
            for ( UInt k = 0; k < nDimensions; ++k )
            {
                for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s1 +=  pow ( Jk[ ig ], (-2.0 / 3.0) ) * Fk[ icoor ][  k ][ ig ] *
                           fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );

                    s2 +=  1.0 / 3.0 * ( Ic_isok[ ig ] * ( 1 / Jk[ig] ) ) *
                           CofFk[ icoor ][ k ][ ig ] * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                }
            }
            vec ( i ) += (s1 - s2) * coef;
        }
    }
}
//! -----------------------------------------------------------------------------------------------






//! Jacobian matrix isochoric part ----------------------------------------------------------------

//! 1. Jacobian matrix : Int { -2/3 * coef * J^(-5/3) *( CofF : \nabla \delta ) ( F : \nabla \v ) }
void stiff_Jac_P1iso_NH_1term ( Real coef,
                                const boost::multi_array<Real, 3 >& CofFk,
                                const boost::multi_array<Real, 3 >& Fk,
                                const std::vector<Real>& Jk ,
                                MatrixElemental& elmat,
                                const CurrentFE& fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt l = 0; l < nDimensions; ++l )
                    {
                        for ( UInt k = 0; k < nDimensions; ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += pow ( Jk[ig], -5. / 3. ) *
                                     Fk[ jcoor ][ l ][ ig ] * fe.phiDer ( j, l, ig ) *
                                     CofFk[ icoor ][ k ][ ig ] * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}



//! 2. Stiffness matrix: Int { 2/9 * coef * ( Ic_iso / J^2 )( CofF : \nabla \delta ) ( CofF : \nabla \v ) }
void stiff_Jac_P1iso_NH_2term ( Real coef,
                                const boost::multi_array<Real, 3 >& CofFk,
                                const std::vector<Real>& Jk ,
                                const std::vector<Real>& Ic_isok,
                                MatrixElemental& elmat,
                                const CurrentFE& fe )
{
    Real s;

    for ( UInt icoor = 0; icoor <  nDimensions; ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt l = 0; l < nDimensions; ++l )
                    {
                        for ( UInt k = 0; k < nDimensions; ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( 1 / (Jk[ig] * Jk[ig] ) ) *  Ic_isok[ig] *
                                     CofFk[ jcoor ][ l ][ ig ] * fe.phiDer ( j, l, ig ) *
                                     CofFk[ icoor ][ k ][ ig ] * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}



//! 3. Stiffness matrix : Int { coef * J^(-2/3) (\nabla \delta : \nabla \v)}
void stiff_Jac_P1iso_NH_3term ( Real          coef,
                                const std::vector<Real>&   Jk,
                                MatrixElemental& elmat,
                                const CurrentFE& fe )
{
    Real s;

    //! assembling diagonal block
    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0.0;
            for ( UInt k = 0; k < nDimensions; ++k )
            {
                for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += pow ( Jk[ig], -2. / 3.) * fe.phiDer ( i, k, ig ) *
                         fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        //! copy of diagonal block
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }
}

//! 4. Stiffness matrix : Int { -2/3 * coef * J^(-5/3) ( F : \nabla \delta ) ( CofF : \nabla \v ) }
void stiff_Jac_P1iso_NH_4term ( Real coef,
                                const boost::multi_array<Real, 3 >& CofFk,
                                const boost::multi_array<Real, 3 >& Fk,
                                const std::vector<Real>& Jk ,
                                MatrixElemental& elmat,
                                const CurrentFE& fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt l = 0; l < nDimensions; ++l )
                    {
                        for ( UInt k = 0; k < nDimensions; ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += pow ( Jk[ig], -5. / 3. ) *
                                     Fk[ icoor ][ k ][ ig ]  * fe.phiDer ( i, k, ig ) *
                                     CofFk[ jcoor ][ l ][ ig ] * fe.phiDer ( j, l, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}



//! 5. Stiffness matrix : Int { 1/3 * coef * J^(-2) * Ic_iso * (CofF [\nabla \delta]^t CofF ) : \nabla \v }
void stiff_Jac_P1iso_NH_5term ( Real coef,
                                const boost::multi_array<Real, 3 >& CofFk,
                                const std::vector<Real>& Jk ,
                                const std::vector<Real>& Ic_isok,
                                MatrixElemental& elmat,
                                const CurrentFE& fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt l = 0; l < nDimensions; ++l )
                    {
                        for ( UInt k = 0; k < nDimensions; ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( 1 / ( Jk[ig] * Jk[ig] ) ) * Ic_isok[ig] *
                                     CofFk[ icoor ][ l ][ ig ] * fe.phiDer ( j, l, ig ) *
                                     CofFk[ jcoor ][ k ][ ig ] * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += s * coef;
                }
            }
        }
    }
}
//! ***********************************************************************************************
//! END OF NEO-HOOKEAN MODEL
//! ***********************************************************************************************

//! ***********************************************************************************************
//! METHODS FOR EXPONENTIAL MODEL
//! ***********************************************************************************************

//! Stiffness vector isochoric part ---------------------------------------------------------------

// Source term : Int { coef * exp(coefExp *(  Ic_iso -3 )) * ( J^(-2/3)* (F : \nabla v) - 1/3 * (Ic_iso / J) * (CofF : \nabla v) ) }
void  source_P1iso_Exp ( Real             coef,
                         Real             coefExp,
                         const boost::multi_array<Real, 3 >& CofFk,
                         const boost::multi_array<Real, 3 >& Fk,
                         const std::vector<Real>&   Jk,
                         const std::vector<Real>&   Ic_isok,
                         VectorElemental& elvec,
                         const CurrentFE& fe )
{

    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        VectorElemental::vector_view vec =  elvec.block( icoor );
        for( UInt i = 0; i < fe.nbFEDof(); ++i )
        {
            s = 0.0;
            for ( UInt k = 0; k < nDimensions; ++k )
            {
                for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += exp ( coefExp * ( Ic_isok[ ig ] - 3.0 ) ) *
                         (pow ( Jk[ ig ], (-2.0 / 3.0) ) * Fk[ icoor ][  k ][ ig ] -
                          1.0 / 3.0 * ( 1 / Jk[ ig ] ) * Ic_isok[ ig ] *
                          CofFk[ icoor ][ k ][ ig ] ) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );

                }
            }
            vec( i ) += s * coef;
        }
    }
}




//! Jacobian matrix isochoric part ----------------------------------------------------------------

//! 1. Stiffness term : Int { - 2/3 *coef * J^(-5/3) * exp( coefExp*( Ic_iso - 3) )* ( 1. + coefExp * Ic_iso ) * ( CofF : \nabla \delta ) ( F : \nabla \v ) }
void  stiff_Jac_P1iso_Exp_1term ( Real             coef,
                                  Real             coefExp,
                                  const boost::multi_array<Real, 3 >& CofFk,
                                  const boost::multi_array<Real, 3 >& Fk,
                                  const std::vector<Real>&   Jk ,
                                  const std::vector<Real>&   Ic_isok,
                                  MatrixElemental&         elmat,
                                  const CurrentFE& fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt l = 0; l < nDimensions; ++l )
                    {
                        for ( UInt k = 0; k < nDimensions; ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s +=std::pow( Jk[ig], -5./3. ) * exp( coefExp*( Ic_isok[ig] - 3 ) ) *
                                    ( 1. + coefExp * Ic_isok[ig] ) *
                                    CofFk[ jcoor ][ k ][ ig ] * fe.phiDer( j, k, ig ) *
                                    Fk[ icoor ][ l ][ ig ] * fe.phiDer( i, l, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat ( i, j ) +=  coef * s;
                }
            }

	    }
	}
}

//! 2. Stiffness term : Int { 2 * coef * coefExp * J^(-4/3) * exp( coefExp*( Ic_iso - 3) ) * ( F : \nabla \delta ) ( F : \nabla \v )}
void  stiff_Jac_P1iso_Exp_2term ( Real             coef,
                                  Real             coefExp,
                                  const boost::multi_array<Real, 3 >& Fk,
                                  const std::vector<Real>&   Jk,
                                  const std::vector<Real>&   Ic_isok,
                                  MatrixElemental&         elmat,
                                  const CurrentFE& fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt l = 0; l < nDimensions; ++l )
                    {
                        for ( UInt k = 0; k < nDimensions; ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s +=std::pow( Jk[ig], -4.0/3.0 ) * exp( coefExp*(  Ic_isok[ig] -3  ) ) *
                                    Fk[ jcoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) *
                                    Fk[ icoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}

//! 3. Stiffness term: Int { 2.0/9.0 * coef * J^-2 * Ic_iso * exp( coefExp*( Ic_iso - 3) ) * ( 1. + coefExp * Ic_iso )( CofF : \nabla \delta ) ( CofF : \nabla \v )}
void  stiff_Jac_P1iso_Exp_3term ( Real coef, Real  coefExp,
                                  const boost::multi_array<Real, 3 >& CofFk,
                                  const std::vector<Real>&   Jk,
                                  const std::vector<Real>&   Ic_isok,
                                  MatrixElemental&         elmat,
                                  const CurrentFE& fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt l = 0; l < nDimensions; ++l )
                    {
                        for ( UInt k = 0; k < nDimensions; ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( 1 / ( Jk[ig] * Jk[ig] ) ) * exp ( coefExp * ( Ic_isok[ig] - 3 ) ) *
                                     ( 1. + coefExp * Ic_isok[ig] ) * Ic_isok[ig] *
                                     CofFk[ jcoor ][ l ][ ig ] * fe.phiDer ( j, l, ig ) *
                                     CofFk[ icoor ][ k ][ ig ] * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );

                            }
                        }
                    }
                    mat ( i, j ) +=  coef * s;
                }
            }
        }
    }
}

//! 4. Stiffness term: Int { -2.0/3.0 * coef * J^(-5/3) * exp( coefExp*( Ic_iso - 3) ) * ( 1. + coefExp * Ic_iso )( F : \nabla \delta ) ( CofF : \nabla \v ) }
void  stiff_Jac_P1iso_Exp_4term ( Real coef, Real  coefExp,
                                  const boost::multi_array<Real, 3 >& CofFk,
                                  const boost::multi_array<Real, 3 >& Fk,
                                  const std::vector<Real>&   Jk ,
                                  const std::vector<Real>&   Ic_isok,
                                  MatrixElemental&         elmat,
                                  const CurrentFE& fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt l = 0; l < nDimensions; ++l )
                    {
                        for ( UInt k = 0; k < nDimensions; ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s +=std::pow( Jk[ig], -5./3. ) * exp( coefExp*( Ic_isok[ig] - 3  ) ) *
                                    ( 1. + coefExp * Ic_isok[ig] ) *
                                    Fk[ jcoor ][ k ][ ig ]  * fe.phiDer( j, k, ig ) *
                                    CofFk[ icoor ][ l ][ ig ] *  fe.phiDer( i, l, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += coef * s;

                }
            }
        }
    }
}

//! 5. Stiffness term : Int {coef * J^(-2/3) * exp( coefExp*( Ic_iso - 3)) (\nabla \delta: \nabla \v)}
void  stiff_Jac_P1iso_Exp_5term ( Real             coef,
                                  Real             coefExp,
                                  const std::vector<Real>&   Jk,
                                  const std::vector<Real>&   Ic_isok,
                                  MatrixElemental&         elmat,
                                  const CurrentFE& fe )
{
    Real s;

    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );
    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0.0;
            for ( UInt k = 0; k < nDimensions; ++k )
            {
                for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s +=std::pow(Jk[ig], -2.0/3.0) * exp( coefExp*( Ic_isok[ig] -3  ) ) *
                        fe.phiDer( i, k, ig ) *  fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        //! copy of diagonal block
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }
}

//! 6. Stiffness term : Int { 1.0/3.0 * coef * J^(-2) * Ic_iso *  exp(coefExp( Ic_iso - 3)) * (CofF [\nabla \delta]^t CofF ) : \nabla \v }
void  stiff_Jac_P1iso_Exp_6term ( Real             coef,
                                  Real             coefExp,
                                  const boost::multi_array<Real, 3 >& CofFk,
                                  const std::vector<Real>&   Jk,
                                  const std::vector<Real>&   Ic_isok,
                                  MatrixElemental&         elmat,
                                  const CurrentFE& fe )
{
    Real s;

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt l = 0; l < nDimensions; ++l )
                    {
                        for ( UInt k = 0; k < nDimensions; ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( 1/( Jk[ig]*Jk[ig] ) ) * Ic_isok[ig] *
                                    exp( coefExp*( Ic_isok[ig] -3  ) ) *
                                    CofFk[ icoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) *
                                    CofFk[ jcoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat ( i, j ) +=  s * coef;
                }
            }
        }
    }
}

//! ***********************************************************************************************
//! END OF EXPONENTIAL MODEL
//! ***********************************************************************************************

//! ***********************************************************************************************
//! METHODS FOR TENSORIAL CALCULUS
//! ***********************************************************************************************

//! Computation of the Right Cauchy Green tensor given the tensor F.
void computeInvariantsRightCauchyGreenTensor(std::vector<LifeV::Real>& invariants,
                                             const Epetra_SerialDenseMatrix& tensorF,
                                             Epetra_SerialDenseMatrix& cofactorF)
{

    //Computation of the invariants
    //At the moment, only the first one is really computed.
    //The others are not still used in the constitutive laws

    Real C11(0);
    Real C22(0);
    Real C33(0);

    //It is not rescaled by the determinant. It is done inside the method to compute the local Piola
    //cofactorF.Scale(invariants[3]);

    C11 = tensorF(0,0)*tensorF(0,0) + tensorF(1,0)*tensorF(1,0) + tensorF(2,0)*tensorF(2,0);
    C22 = tensorF(0,1)*tensorF(0,1) + tensorF(1,1)*tensorF(1,1) + tensorF(2,1)*tensorF(2,1);
    C33 = tensorF(0,2)*tensorF(0,2) + tensorF(1,2)*tensorF(1,2) + tensorF(2,2)*tensorF(2,2);

    invariants[0]=C11 + C22 + C33; //First invariant C
    invariants[1]=0.0; //Second invariant C
    invariants[2]=0.0; //Third invariant C
    invariants[3]=tensorF(0,0) * ( tensorF(1,1)*tensorF(2,2) - tensorF(1,2)*tensorF(2,1) ) - tensorF(0,1) * ( tensorF(1,0)*tensorF(2,2) - tensorF(1,2)*tensorF(2,0) ) + tensorF(0,2) * ( tensorF(1,0)*tensorF(2,1) - tensorF(1,1)*tensorF(2,0) ); //Determinant F

    //Computation of the Cofactor of F
    cofactorF( 0 , 0 ) =   ( tensorF(1,1)*tensorF(2,2) - tensorF(1,2)*tensorF(2,1) );
    cofactorF( 0 , 1 ) = - ( tensorF(1,0)*tensorF(2,2) - tensorF(2,0)*tensorF(1,2) );
    cofactorF( 0 , 2 ) =   ( tensorF(1,0)*tensorF(2,1) - tensorF(1,1)*tensorF(2,0) );
    cofactorF( 1 , 0 ) = - ( tensorF(0,1)*tensorF(2,2) - tensorF(0,2)*tensorF(2,1) );
    cofactorF( 1 , 1 ) =   ( tensorF(0,0)*tensorF(2,2) - tensorF(0,2)*tensorF(2,0) );
    cofactorF( 1 , 2 ) = - ( tensorF(0,0)*tensorF(2,1) - tensorF(2,0)*tensorF(0,1) );
    cofactorF( 2 , 0 ) =   ( tensorF(0,1)*tensorF(1,2) - tensorF(0,2)*tensorF(1,1) );
    cofactorF( 2 , 1 ) = - ( tensorF(0,0)*tensorF(1,2) - tensorF(0,2)*tensorF(1,0) );
    cofactorF( 2 , 2 ) =   ( tensorF(0,0)*tensorF(1,1) - tensorF(1,0)*tensorF(0,1) );

    cofactorF.Scale(1/invariants[3]);
}

//! Computation of the Right Cauchy Green tensor given the tensor F.
void computeInvariantsRightCauchyGreenTensor(std::vector<LifeV::Real>& invariants,
                                             const Epetra_SerialDenseMatrix& tensorF)
{

    invariants[0]=0.0; //First invariant C
    invariants[1]=0.0; //Second invariant C
    invariants[2]=0.0; //Third invariant C
    invariants[3]=tensorF(0,0) * ( tensorF(1,1)*tensorF(2,2) - tensorF(1,2)*tensorF(2,1) ) - tensorF(0,1) * ( tensorF(1,0)*tensorF(2,2) - tensorF(1,2)*tensorF(2,0) ) + tensorF(0,2) * ( tensorF(1,0)*tensorF(2,1) - tensorF(1,1)*tensorF(2,0) ); //Determinant F

}


void computeCauchyStressTensor(Epetra_SerialDenseMatrix& cauchy,
                               Epetra_SerialDenseMatrix& firstPiola,
                               LifeV::Real det,
                               Epetra_SerialDenseMatrix& tensorF)
{

    firstPiola.Scale( 1/det );
    cauchy.Multiply('N','T',1.0,firstPiola,tensorF,0.0);

}

void computeEigenvalues(const Epetra_SerialDenseMatrix& cauchy,
                        std::vector<LifeV::Real>& eigenvaluesR,
                        std::vector<LifeV::Real>& eigenvaluesI)

{

    // LAPACK wrapper of Epetra
    Epetra_LAPACK lapack;

    //List of flags for Lapack Function
    //For documentation, have a look at http://www.netlib.org/lapack/double/dgeev.f

    char JOBVL = 'N';
    char JOBVR = 'N';

    //Size of the matrix
    Int Dim = cauchy.RowDim();

    //Arrays to store eigenvalues (their number = nDimensions)
    double WR[nDimensions];
    double WI[nDimensions];

    //Number of eigenvectors
    Int LDVR = nDimensions;
    Int LDVL = nDimensions;

    //Arrays to store eigenvectors
    Int length = nDimensions * 3;

    double VR[length];
    double VL[length];

    Int LWORK = 9;
    Int INFO = 0;

    double WORK[LWORK];

    double A[length];

    for (UInt i(0); i< nDimensions; i++)
        for (UInt j(0);j<nDimensions; j++)
            A[nDimensions * i + j] = cauchy(i,j);

    lapack.GEEV(JOBVL, JOBVR, Dim, A /*cauchy*/, Dim, &WR[0], &WI[0], VL, LDVL, VR, LDVR, WORK, LWORK, &INFO);
    ASSERT_PRE( !INFO, "Calculation of the Eigenvalues failed!!!" );

    for( UInt i(0); i < nDimensions; i++ )
    {
        eigenvaluesR[i] = WR[i];
        eigenvaluesI[i] = WI[i];
    }

}

//! ***********************************************************************************************
//! METHODS FOR THE ST. VENANT KIRCHHOFF PENALIZED LAW
//! ***********************************************************************************************
//! Stiffness vector isochoric part ---------------------------------------------------------------

// Source term : Int { ( \frac{lambda}{2} * Ic_iso - \frac{3}{2}*lambda - mu ) (F : \nabla v) - 1/3 * (Ic) * (F^-T : \nabla v) ) }
void  source_P1iso_VKPenalized( Real             lambda,
                                Real             mu,
                                const boost::multi_array<Real,3 >& FkMinusTransposed,
                                const boost::multi_array<Real,3 >& Fk,
                                const std::vector<Real>&   Ic_isok,
                                const std::vector<Real>&   Ic_k,
                                const std::vector<Real>&   Jack_k,
                                VectorElemental& elvec,
                                const CurrentFE& fe )
{

    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        VectorElemental::vector_view vec =  elvec.block( icoor );
        for( UInt i = 0; i < fe.nbFEDof(); ++i )
        {
            s = 0.0;
            for( UInt k = 0; k < nDimensions; ++k )
            {
                for( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s +=  std::pow( Jack_k[ ig ], (-2.0/3.0) )  * ( ( lambda/2.0 ) * Ic_isok[ ig ] - (3.0/2.0)* lambda - mu  ) *
                        (Fk[ icoor ][  k ][ ig ] - (1.0/3.0) * Ic_k[ ig ] * FkMinusTransposed[ icoor ][ k ][ ig ] ) * fe.phiDer( i, k, ig ) * fe.weightDet( ig );

                }
            }
            vec( i ) += s;
        }
    }
}


// Source term : Int { ( 2 * mu * Jk^(-4.0/3.0) ) * ( (F*C : \nabla v) - 1/3 * (Ic_Squared) * (F^-T : \nabla v) ) }
void  source_P2iso_VKPenalized( Real             mu,
                                const boost::multi_array<Real,3 >& FkMinusTransposed,
                                const boost::multi_array<Real,3 >& FkCk,
                                const std::vector<Real>&   Ic_Squared,
                                const std::vector<Real>&   Jk,
                                VectorElemental& elvec,
                                const CurrentFE& fe )
{

    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        VectorElemental::vector_view vec =  elvec.block( icoor );
        for( UInt i = 0; i < fe.nbFEDof(); ++i )
        {
            s = 0.0;
            for( UInt k = 0; k < nDimensions; ++k )
            {
                for( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += ( mu * std::pow( Jk[ ig ], (-4.0/3.0) ) ) * (FkCk[ icoor ][  k ][ ig ] - (1.0/3.0)  * Ic_Squared[ ig ] * FkMinusTransposed[ icoor ][ k ][ ig ] )* fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                }
            }
            vec( i ) += s;
        }
    }
}

//! Jacobian matrix of the first Piola-Kirchhoff tensor for the VK Penalized law
//! 1. Stiffness term : int { -(2.0/3.0) * Jk^(-2.0/3.0) * ( (lambda/2) * Ic_isok - ( (3/2)*lambda + mu ) ) * F^-T:\nabla \delta ) * ( F - (1.0/3.0) * Ic_k * F^-T ): \nabla \v  }
void  stiff_Jac_P1iso_VKPenalized_0term( Real             lambda, Real mu,
                                         const boost::multi_array<Real,3 >& FkMinusTransposed,
                                         const boost::multi_array<Real,3 >& Fk,
                                         const std::vector<Real>&   Jk ,
                                         const std::vector<Real>&   Ic_k ,
                                         const std::vector<Real>&   IcIso_k ,
                                         MatrixElemental&         elmat,
                                         const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += (-2.0/3.0) * std::pow(Jk[ ig ], (-2.0/3.0) ) * ( ( (lambda/2.0) * IcIso_k[ ig ] ) - ( (3.0/2.0) * lambda + mu)  ) *
                                    fe.phiDer( i, l, ig ) * ( Fk[ icoor ][ l ][ ig ] - (1.0/3.0) * Ic_k[ ig ] * FkMinusTransposed[ icoor ][ l ][ ig ] ) *
                                    FkMinusTransposed[ jcoor ][ k ][ ig ] * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) +=  s;
                }
            }
	    }
	}
}

//! 1. Stiffness term :int{ J^(-2/3) * (lambda / 2) * ( (-2/3) * Ic_k * J^(-2/3) * F^-T:\nabla \delta ) * ( F : \nabla \v )}
void  stiff_Jac_P1iso_VKPenalized_1term( Real             coeff,
                                         const boost::multi_array<Real,3 >& FkMinusTransposed,
                                         const boost::multi_array<Real,3 >& Fk,
                                         const std::vector<Real>&   Jk ,
                                         const std::vector<Real>&   Ic_k ,
                                         MatrixElemental&         elmat,
                                         const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( -2.0/3.0 ) * Ic_k[ ig ] *std::pow( Jk[ig], (-4.0/3.0) ) *
                                    fe.phiDer( i, l, ig ) * Fk[ icoor ][ l ][ ig ] *
                                    FkMinusTransposed[ jcoor ][ k ][ ig ] * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) +=  coeff * s;
                }
            }
	    }
	}
}

//! 3. Stiffness term:int { J^(-2/3) * (lambda / 2) * ( ( 2/9 ) * J^(-2/3) * Ic_k^2 ) * ( F^-T : \nabla \delta ) ( F^-T : \nabla \v ) }
void  stiff_Jac_P1iso_VKPenalized_2term( Real coef,
                                         const boost::multi_array<Real,3 >& FkMinusTransposed,
                                         const std::vector<Real>&   Jk,
                                         const std::vector<Real>&   Ic_k,
                                         MatrixElemental&         elmat,
                                         const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( ( 2.0 / 9.0 ) * std::pow( Jk[ig], (-4.0/3.0) ) * Ic_k[ig] * Ic_k[ig] ) *
                                    FkMinusTransposed[ icoor ][ l ][ ig ] * fe.phiDer( i, l, ig ) *
                                    FkMinusTransposed[ jcoor ][ k ][ ig ] * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) +=  coef * s;
                }
            }
	    }
	}
}

//! 2. Stiffness term : Int { J^(-2/3) * coeff * ( 2 * J^(-2/3) ) * ( F : \nabla \delta ) ( F : \nabla \v )}
void  stiff_Jac_P1iso_VKPenalized_3term( Real             coef,
                                         const boost::multi_array<Real,3 >& Fk,
                                         const std::vector<Real>&   Jk,
                                         MatrixElemental&         elmat,
                                         const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s +=  2.0 *std::pow( Jk[ig], (-4.0/3.0) ) *
                                    fe.phiDer( i, l, ig ) * Fk[ icoor ][ l ][ ig ] *
                                    Fk[ jcoor ][ k ][ ig ] * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}

//! 4. Stiffness term:int{J^(-2/3) * (lambda/2) * ( -2.0/3.0 * J^(-2/3) * Ic_k ) * ( F : \nabla \delta)*(F^-T : \nabla \v )}
void  stiff_Jac_P1iso_VKPenalized_4term( Real coef,
                                         const boost::multi_array<Real,3 >& FkMinusTransposed,
                                         const boost::multi_array<Real,3 >& Fk,
                                         const std::vector<Real>&   Jk ,
                                         const std::vector<Real>&   Ic_k,
                                         MatrixElemental&         elmat,
                                         const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( ( -2.0/3.0 ) * std::pow( Jk[ig], (-4.0/3.0) ) * Ic_k[ig] ) *
                                    fe.phiDer( i, l, ig ) * FkMinusTransposed[ icoor ][ l ][ ig ] *
                                    Fk[ jcoor ][ k ][ ig ]  * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}

//! 5. Stiffness term : int {  J^(-2.0/3.0) * ( (lambda/2) * Ic_isok - ( (3/2)*lambda + mu ) ) * \nabla \delta : \nabla v }
void  stiff_Jac_P1iso_VKPenalized_5term( Real             coef,
                                         Real             secondCoef,
                                         const std::vector<Real>&   Jk,
                                         const std::vector<Real>&   Ic_isok,
                                         MatrixElemental&         elmat,
                                         const CurrentFE& fe )
{
    Real s;

    MatrixElemental::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    for( UInt i = 0; i < fe.nbFEDof(); ++i )
	{
        for( UInt j = 0; j < fe.nbFEDof(); ++j )
	    {
            s = 0.0;
            for( UInt k = 0; k < nDimensions; ++k )
            {
                for( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += std::pow( Jk[ ig ], (-2.0/3.0) ) * ( (coef/2.0) * Ic_isok[ ig ] - ( (3.0/2.0) * coef + secondCoef ) ) *
                        fe.phiDer( i, k, ig ) *  fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                }
	    	}
            mat_tmp( i, j ) = s;
	    }
	}

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        //! copy of diagonal block
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
	}
}

//! 6. Stiffness term: int { J^(-2.0/3.0) * ( (lambda/2) * Ic_isok - ( (3/2)*lambda + mu ) ) * ( (-2/3) * ( F :\nabla \delta ) ) * ( F^-T : \nabla v ) }
void  stiff_Jac_P1iso_VKPenalized_6term( Real coef, Real  secondCoef,
                                         const std::vector<Real>&   Jk,
                                         const std::vector<Real>&   Ic_isok,
                                         const boost::multi_array<Real,3 >& Fk,
                                         const boost::multi_array<Real,3 >& FkMinusTransposed,
                                         MatrixElemental&         elmat,
                                         const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += std::pow( Jk[ ig ], (-2.0/3.0) ) * (-2.0/3.0) * ( (coef/2.0) * Ic_isok[ ig ] - ( (3.0/2.0) * coef + secondCoef ) ) *
                                    fe.phiDer( i, l, ig ) * FkMinusTransposed[ icoor ][ l ][ ig ] *
                                    Fk[ jcoor ][ k ][ ig ]  * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += s;
                }
            }
	    }
	}
}

//! 7. Stiffness term : int { ( J^(-2/3) * (lambda/2) * Ic_isok - ( (3/2)*lambda + mu ) ) * ( (1/3) * Ic_k * ( F^-T \nabla \delta^T F-T ) : \nabla v  }
void  stiff_Jac_P1iso_VKPenalized_7term( Real             coef,
                                         Real             secondCoef,
                                         const boost::multi_array<Real,3 >& FkMinusTransposed,
                                         const std::vector<Real>&   Ic_isok,
                                         const std::vector<Real>&   Ic_k,
                                         const std::vector<Real>&   Jk,
                                         MatrixElemental&         elmat,
                                         const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += std::pow( Jk[ ig ], (-2.0/3.0) ) * ( (coef/2.0) * Ic_isok[ ig ] - ( (3.0/2.0) * coef + secondCoef ) ) * ( (1.0/3.0) * Ic_k[ ig ] ) *
                                    FkMinusTransposed[ icoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) *
                                    FkMinusTransposed[ jcoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) +=  s;
                }
            }
	    }
	}
}

//! 8. Stiffness term : int { ( -4.0/3.0) * ( mu * J^(-4/3) ) * ( F^-T: \grad \delta ) * ( F C ) : \nabla v  }
void  stiff_Jac_P1iso_VKPenalized_8term( Real coef, const std::vector<Real>& Jack_k,
                                         const boost::multi_array<Real,3 >& FkMinusTransposed,
                                         const boost::multi_array<Real,3 >& FkCk,
                                         MatrixElemental& elmat,
                                         const CurrentFE& fe )
{

    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += (-4.0/3.0) *  std::pow( Jack_k[ ig ], (-4.0/3.0) ) *
                                    fe.phiDer( i, l, ig ) * FkCk[ icoor ][ l ][ ig ] *
                                    FkMinusTransposed[ jcoor ][ k ][ ig ]  * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += s * coef;
                }
            }
	    }
	}
}


//! 8. Stiffness term : int { ( 4.0/9.0) * ( mu * J^(-4/3) ) Ic_kSquared * (F^-T : \grad \delta ) * F^-T : \nabla \v  }
void  stiff_Jac_P1iso_VKPenalized_9term( Real coef, const std::vector<Real>& Jack_k,
                                         const std::vector<Real>& Ic_kSquared,
                                         const boost::multi_array<Real,3 >& FkMinusTransposed,
                                         MatrixElemental& elmat,
                                         const CurrentFE& fe )
{
    Real s;
    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( (-4.0/3.0) * std::pow( Jack_k[ ig ], (-4.0/3.0) ) ) *
                                    fe.phiDer( i, l, ig ) * (-1.0/3.0) * Ic_kSquared[ ig ] * FkMinusTransposed[ icoor ][ l ][ ig ] *
                                    FkMinusTransposed[ jcoor ][ k ][ ig ]  * fe.phiDer( j, k, ig ) *fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += s * coef;
                }
            }
	    }
	}
}




//! 10. Stiffness term : int { ( mu * J^(-4/3) ) * ( \nabla \delta * C ) : \nabla v  }
void  stiff_Jac_P1iso_VKPenalized_10term( Real             coef,
                                         const std::vector<Real>&   Jack_k,
                                         const boost::multi_array<Real,3 >& Ck,
                                         MatrixElemental&         elmat,
                                         const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
        for( UInt i = 0; i < fe.nbFEDof(); ++i )
        {
            for( UInt j = 0; j < fe.nbFEDof(); ++j )
            {
                s = 0.0;
                for( UInt l = 0; l < nDimensions; ++l )
                {
                    for( UInt k = 0; k < nDimensions; ++k )
                    {
                        for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                        {
                            s += std::pow( Jack_k[ ig ], (-4.0/3.0) ) * Ck[ l ][ k ][ ig ] * fe.phiDer( i, k, ig ) *
                                fe.phiDer( j, l, ig ) * fe.weightDet( ig );
                        }
                    }
                }
                mat( i, j ) +=  coef * s;
            }
        }
    }
}

//! 6. Stiffness term : int { ( mu * J^(-4/3) ) * (F [\nabla \delta]^T F ) : \nabla \v  }
void  stiff_Jac_P1iso_VKPenalized_11term( Real             coef,
                                         const std::vector<Real>&   Jk,
                                         const boost::multi_array<Real,3 >& Fk,
                                         MatrixElemental&         elmat,
                                         const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += std::pow( Jk[ ig ], (-4.0/3.0)) * fe.phiDer( i, l, ig ) * Fk[ l ][ jcoor ][ ig ] *
                                    Fk[ icoor ][ k ][ ig ] *  fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) +=  coef * s;
                }
            }
	    }
	}
}

//! 6. Stiffness term : int { ( mu * J^(-4/3) ) * (F * F^T * [\nabla \delta] ) : \nabla \v  }
void  stiff_Jac_P1iso_VKPenalized_12term( Real             coef,
                                          const std::vector<Real>&   Jk,
                                          const boost::multi_array<Real,3 >& Fk,
                                          MatrixElemental&         elmat,
                                          const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt p = 0; p < nDimensions; p++ )
                    {
                        for( UInt k = 0; k < nDimensions; k++ )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += std::pow( Jk[ig], (-4.0/3.0) ) * Fk[ icoor ][ p ][ ig ] * fe.phiDer( i, k, ig ) *
                                    Fk[ jcoor ][ p ][ ig ] * fe.phiDer( j, k, ig ) *  fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) +=  coef * s;
                }
            }
	    }
	}
}

//! 11. Stiffness term : int {  ( mu * J^(-4/3) * ( (1/3) *  Ic_SquaredK * ( F^-T [\nabla \delta ]^T F^-T) : \nabla v ) }
void  stiff_Jac_P1iso_VKPenalized_13term( Real             coef,
                                          const std::vector<Real>&   Jk,
                                          const std::vector<Real>&   Ic_kSquared,
                                          const boost::multi_array<Real,3 >& FkMinusTransposed,
                                          MatrixElemental&         elmat,
                                          const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( std::pow( Jk[ig], (-4.0/3.0) ) ) * ( Ic_kSquared[ig] / 3.0 ) *
                                    fe.phiDer( j, l, ig ) * FkMinusTransposed[ icoor ][ l ][ ig ] *
                                    fe.phiDer( i, k, ig ) * FkMinusTransposed[ jcoor ][ k ][ ig ] * fe.weightDet( ig );\
                            }
                        }
                    }
                    mat( i, j ) +=  s * coef;
                }
            }
	    }
	}
}

//! 12. Stiffness term : int {  ( mu * J^(-4/3) ) * ( (-4.0/3.0) * ( FkCk : \nabla \delta ) ) * F^-T : \nabla v ) }
void  stiff_Jac_P1iso_VKPenalized_14term( Real             coef,
                                          const std::vector<Real>&   Jk,
                                          const boost::multi_array<Real,3 >& FkCk,
                                          const boost::multi_array<Real,3 >& FkMinusTransposed,
                                          MatrixElemental&         elmat,
                                          const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( std::pow( Jk[ig], (-4.0/3.0) ) ) * ( -4.0 / 3.0 ) *
                                    fe.phiDer( i, l, ig ) * FkMinusTransposed[ icoor ][ l ][ ig ] *
                                    ( FkCk[ jcoor ][ k ][ ig ] * fe.phiDer( j, k, ig ) ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) +=  s * coef;
                }
            }
	    }
	}
}

//! ***********************************************************************************************
//! END OF VENANT-KIRCHHOFF PENALIZED MODEL
//! ***********************************************************************************************

//! ***********************************************************************************************
//! SECOND ORDER EXPONENTIAL MODEL
//! ***********************************************************************************************
//! Methods for the Second Order Exponential law
//! Assemble the first Piola-Kirchhoff tensor
//int { 2 * alpha * ( Ic1_iso - 3 ) * exp(gamma *(  Ic1_iso -3 )^2) * ( J1^(-2/3)* (F1 : \nabla v) - 1/3 * (Ic1_iso / J1) * (CofF1 : \nabla v) ) }
void  source_P1iso_SecondOrderExponential( Real coef, Real coefExp,
                                           const boost::multi_array<Real,3 >& CofFk,
                                           const boost::multi_array<Real,3 >& Fk,
                                           const std::vector<Real>&   Jk,
                                           const std::vector<Real>&   Ic_isok,
                                           VectorElemental& elvec,
                                           const CurrentFE& fe )
{

    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
    {
        VectorElemental::vector_view vec =  elvec.block( icoor );
        for( UInt i = 0; i < fe.nbFEDof(); ++i )
        {
            s = 0.0;
            for( UInt k = 0; k < nDimensions; ++k )
            {
                for( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += (Ic_isok[ ig ] - 3.0) * exp( coefExp * (Ic_isok[ ig ] - 3.0) * (Ic_isok[ ig ] - 3.0) ) *
                        (std::pow( Jk[ ig ], (-2.0/3.0) ) * Fk[ icoor ][  k ][ ig ]
                         - 1.0/3.0 * ( Ic_isok[ ig ] / Jk[ ig ] ) * CofFk[ icoor ][ k ][ ig ] ) * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                }
            }
            vec( i ) += s * coef;
        }
    }
}

//! 1. Stiffness term : Int { ( - 4/3 *coef * J^(-5/3) * exp( coefExp*( Ic_iso - 3)^2 ) * ( tr_isoCk + ( tr_isoCk - 3 ) * ( 2 *coefExp * ( tr_isoCk - 3 ) + 1 ) ) ) * ( CofF : \nabla \delta ) ( F : \nabla \v ) }
// In this method, the coef that is passed is -4/3 * alpha, therefore it can multiply the sum at the end.
void  stiff_Jac_P1iso_SecondOrderExp_1term( Real             coef,
                                            Real             coefExp,
                                            const boost::multi_array<Real,3 >& CofFk,
                                            const boost::multi_array<Real,3 >& Fk,
                                            const std::vector<Real>&   Jk ,
                                            const std::vector<Real>&   Ic_isok,
                                            MatrixElemental&         elmat,
                                            const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( std::pow( Jk[ig], -5.0/3.0) * exp(coefExp * (Ic_isok[ig] - 3.0) * (Ic_isok[ig] - 3.0) ) ) *
                                    ( (Ic_isok[ig] - 3.0) * ( 1.0 + 2 * coefExp * (Ic_isok[ig] - 3.0) * Ic_isok[ig] ) + Ic_isok[ig] ) *
                                    fe.phiDer( i, l, ig ) * Fk[ icoor ][ l ][ ig ] *
                                    CofFk[ jcoor ][ k ][ ig ] * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) +=  coef * s;
                }
            }
	    }
	}
}

//! 2. Stiffness term : Int { 4 * coef * J^(-4/3) * exp( coefExp*( Ic_iso - 3)*( Ic_iso - 3) ) *
//!                          ( 1.0 + 2 * coefExp * (Ic_isoK - 3)^2 ) * ( F : \nabla \delta ) ( F : \nabla \v )}
// When the method is called, the coef parameter stores already 4 * alpha. This is why it is used at the end.
void  stiff_Jac_P1iso_SecondOrderExp_2term( Real             coef,
                                            Real             coefExp,
                                            const boost::multi_array<Real,3 >& Fk,
                                            const std::vector<Real>&   Jk,
                                            const std::vector<Real>&   Ic_isok,
                                            MatrixElemental&         elmat,
                                            const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += std::pow( Jk[ig], -4.0/3.0 ) * exp( coefExp * (Ic_isok[ig] - 3.0)*(Ic_isok[ig] - 3.0) ) *
                                    ( 1.0 + 2 * coefExp * (Ic_isok[ig] - 3.0)*(Ic_isok[ig] - 3.0) ) *
                                    fe.phiDer( i, l, ig ) * Fk[ icoor ][ l ][ ig ] *
                                    Fk[ jcoor ][ k ][ ig ] * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}

//! 3. Stiffness term:
//Int { ( 4.0/9.0 *  alpha * J^-2 *  exp( gamma*( Ic_iso - 3)^2 ) * ( Ic_isoK + 2 * gamma * (Ic_isoK - 3)^2 * Ic_isoK + (Ic_isoK - 3) ) ) *
//      ( CofF : \nabla \delta ) ( CofF : \nabla \v )}
// The coef that is passed stores 4/9 * alpha.
void  stiff_Jac_P1iso_SecondOrderExp_3term( Real coef, Real  coefExp,
                                            const boost::multi_array<Real,3 >& CofFk,
                                            const std::vector<Real>&   Jk,
                                            const std::vector<Real>&   Ic_isok,
                                            MatrixElemental&         elmat,
                                            const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s +=( ( 1/(Jk[ig]*Jk[ig]) )* Ic_isok[ig] * exp( coefExp * (Ic_isok[ig] - 3.0) * (Ic_isok[ig] - 3.0) ) ) *
                                    ( (Ic_isok[ig] - 3.0) * ( 1.0 + 2.0 * coefExp * (Ic_isok[ig] - 3.0) * Ic_isok[ig] ) + Ic_isok[ig] ) *
                                    CofFk[ icoor ][ l ][ ig ] * fe.phiDer( i, l, ig ) *
                                    CofFk[ jcoor ][ k ][ ig ] * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) +=  coef * s;
                }
            }
	    }
	}
}


//! 4. Stiffness term:
//Int { (-4.0/3.0 *  alpha * J^(-5/3) * exp( gamma*( Ic_iso - 3)*( Ic_iso - 3) ) * ( Ic_isoK + ( Ic_isok - 3.0 ) * (2*gamma*(Ic_isok - 3)Ic + 1) ) *
//      ( F : \nabla \delta ) ( CofF : \nabla \v ) }
// As in other methods, the coef that is passed is -4.0/3.0 * alpha.
void  stiff_Jac_P1iso_SecondOrderExp_4term( Real coef, Real  coefExp,
                                            const boost::multi_array<Real,3 >& CofFk,
                                            const boost::multi_array<Real,3 >& Fk,
                                            const std::vector<Real>&   Jk ,
                                            const std::vector<Real>&   Ic_isok,
                                            const std::vector<Real>&   Ic_k,
                                            MatrixElemental&         elmat,
                                            const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += std::pow( Jk[ig], (-5.0/3.0) ) * exp( coefExp * (Ic_isok[ig] - 3.0) * (Ic_isok[ig] - 3.0) ) *
                                    ( Ic_isok[ig] + (Ic_isok[ig] - 3.0) * ( 1.0 + 2.0 * coefExp * (Ic_isok[ig] - 3.0) * Ic_isok[ig] ) ) *
                                    fe.phiDer( i, l, ig ) * CofFk[ icoor ][ l ][ ig ] *
                                    Fk[ jcoor ][ k ][ ig ]  * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}

//! 5. Stiffness term :
//Int {coef * J^(-2/3) * exp( coefExp*( Ic_iso - 3)*( Ic_iso - 3)) * ( Ic_iso - 3)* (\nabla \delta: \nabla \v)}
void  stiff_Jac_P1iso_SecondOrderExp_5term( Real             coef,
                                            Real             coefExp,
                                            const std::vector<Real>&   Jk,
                                            const std::vector<Real>&   Ic_isok,
                                            MatrixElemental&         elmat,
                                            const CurrentFE& fe )
{
    Real s;

    MatrixElemental::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    for( UInt i = 0; i < fe.nbFEDof(); ++i )
	{
        for( UInt j = 0; j < fe.nbFEDof(); ++j )
	    {
            s = 0.0;
            for( UInt k = 0; k < nDimensions; ++k )
            {
                for( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += (Ic_isok[ ig ] - 3.0) * exp( coefExp * (Ic_isok[ig] - 3.0) * (Ic_isok[ig] - 3.0) ) *
                        std::pow( Jk[ig], -2.0/3.0 ) * fe.phiDer( i, k, ig ) *  fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                }
	    	}
            mat_tmp( i, j ) = coef * s;
	    }
	}

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        //! copy of diagonal block
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
	}
}

//! 6. Stiffness term : Int { 2.0/3.0 * coef * J^(-2) * Ic_iso *  exp(coefExp( Ic_iso - 3)) * (CofF [\nabla \delta]^t CofF ) : \nabla \v }
void  stiff_Jac_P1iso_SecondOrderExp_6term( Real             coef,
                                            Real             coefExp,
                                            const boost::multi_array<Real,3 >& CofFk,
                                            const std::vector<Real>&   Jk,
                                            const std::vector<Real>&   Ic_isok,
                                            MatrixElemental&         elmat,
                                            const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for( UInt l = 0; l < nDimensions; ++l )
                    {
                        for( UInt k = 0; k < nDimensions; ++k )
                        {
                            for( UInt ig = 0;ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( 1/(Jk[ig]*Jk[ig]) ) * (Ic_isok[ig] - 3.0) *
                                    exp( coefExp * (Ic_isok[ig] - 3.0) * (Ic_isok[ig] - 3.0) ) * Ic_isok[ig]*
                                    fe.phiDer( j, l, ig ) * CofFk[ icoor ][ l ][ ig ] *
                                    fe.phiDer( i, k, ig ) * CofFk[ jcoor ][ k ][ ig ] *  fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) +=  s * coef;
                }
            }
	    }
	}
}

//! ***********************************************************************************************
//! END OF SECOND ORDER EXPONENTIAL MODEL
//! ***********************************************************************************************

} //! End namespace AssemblyElementalStructure

} //! End namespace LifeV


#endif
