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

void computeGradientLocalDisplacement(boost::multi_array<Real, 3>& gradientLocalDisplacement,
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
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
                }
                gradientLocalDisplacement[ icoor ][ jcoor ][ ig ] = s;
            }
	    }
	}
}


// The methods for linear elastic model (stiff_strain and stiff_div) are implemented in AssemblyElemental.cpp

//! Methods for St. Venant Kirchhoff model
//! Methods for the stiffness matrix

//! \f$ coef \cdot ( trace { [\nabla u^k]^T \nabla u }, \nabla\cdot  v  ) \f$
void stiff_derdiv( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
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

            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
                        {
                            s += fe.phiDer( i, icoor, ig ) * gradientLocalDisplacement[ jcoor ][ k ][ ig ]
                                * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                        }
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}



//! \f$ coef \cdot ( [\nabla u^k]^T \nabla u : \nabla v  )\f$
void stiff_dergradbis( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
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

            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            s += fe.phiDer( i, k, ig ) * gradientLocalDisplacement[ jcoor][ icoor ][ ig ] *
                                fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}



//! \f$  coef * ( (\div u_k) \grad u : \grad v  )
void stiff_divgrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;
    // local vector for \div u^k at each quadrature poInt
    //Real duk[ fe.nbQuadPt() ];
    std::vector<Real > duk(fe.nbQuadPt(),0.0);

    // loop on quadrature poInts
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
	{
        s=0;
        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
	    {
            for ( UInt i = 0; i < fe.nbFEDof(); i++ )
            {
                // construction of \div u^k at a quadrature poInt
                s += fe.phiDer( i, icoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
	    }
        duk[ ig ] = s;
	}

    MatrixElemental::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
	{
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
	    {
            s = 0.0;
            for ( UInt k = 0; k < fe.nbCoor(); ++k )
            {
                for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    s += duk[ ig ] * fe.phiDer( i, k, ig ) *  fe.phiDer( j, k, ig ) * fe.weightDet( ig );
            }
            mat_tmp( i, j ) = coef * s;
	    }
	}

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
	{
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
	}

}

//! \f$ coef * ( \grad u_k : \grad u_k) * ( \grad u : \grad v  )
void stiff_gradgrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s,s1;
    //    (\grad u_k : \grad u_k) at each quadrature poInt
    //Real gguk[ fe.nbQuadPt() ];
    std::vector<Real > gguk(fe.nbQuadPt(),0.0);

    // loop on quadrature poInts
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
	{
        s=0;
        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
	    {
            for ( UInt l = 0; l < fe.nbCoor(); l++ )
            {
                s1=0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    s1+= fe.phiDer( i, l, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
                s += s1*s1;
            }
	    }
        gguk[ ig ] = s;
	}

    MatrixElemental::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
	{
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
	    {
            s = 0.0;
            for ( UInt k = 0; k < fe.nbCoor(); ++k )
            {
                for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    s += gguk[ ig ] * fe.phiDer( i, k, ig ) *  fe.phiDer( j, k, ig ) * fe.weightDet( ig );
            }
            mat_tmp( i, j ) = coef * s;
	    }
	}

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
	{
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
	}
}

void stiff_dergrad_gradbis( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
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

            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                    {
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            s += gradientLocalDisplacement[ icoor ][ jcoor ][ ig ] * fe.phiDer( i, k, ig ) *
                                fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                    }
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}



// coef * ( \grad u^k [\grad u]^T : \grad v )
void stiff_dergrad_gradbis_Tr( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
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

            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                    {
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            s += gradientLocalDisplacement[ icoor ][ k ][ ig ]  * fe.phiDer( j, k, ig ) *
                                fe.phiDer( i, jcoor, ig ) * fe.weightDet( ig );
                    }
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}



// coef * ( \grad u^k [\grad u^k]^T \grad u : \grad v )
void stiff_gradgradTr_gradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;

    //! \grad u^k  [\grad u^k]^T  at each quadrature poInt
    //Real guk_gukT[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];
    boost::multi_array<Real, 3> guk_gukT(boost::extents[fe.nbCoor()][fe.nbCoor()][fe.nbQuadPt()]);

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
                            s  += fe.phiDer( i, n, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ] *
                                fe.phiDer( j, n, ig ) * uk_loc.vec() [ j + jcoor * fe.nbFEDof() ];
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

            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            s += fe.phiDer( i, k, ig ) * guk_gukT[ icoor ][ jcoor ][ ig ] *
                                fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}
// End of methods for the stiffness matrix (St. Venant-Kirchhoff material)

//! Methods for the jacobian (St. Venant-Kirchhoff material)

//! \f$ coef \cdot ( [\nabla u]^T \nabla u^k + [\nabla u^k]^T \nabla u : \nabla v  )\f$
void stiff_dergrad( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement,
                    MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
	{
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
	    {

            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                        {
                            s += fe.phiDer( i, k, ig ) *
                                ( gradientLocalDisplacement[ jcoor ][ k ][ ig ] * fe.phiDer( j, icoor, ig )
                                  + gradientLocalDisplacement[ jcoor ][ icoor ][ ig ] * fe.phiDer( j, k, ig ) ) * fe.weightDet( ig );
                        }
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}



// coef * ( (\div u) \grad u_k : \grad v  )
void stiff_divgrad_2( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement,
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
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            s += fe.phiDer( j, jcoor, ig ) *
                                gradientLocalDisplacement[ icoor ][ k ][ ig ] *
                                fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}

// coef * ( \grad u_k : \grad u) *( \grad u_k : \grad v  )
void stiff_gradgrad_2( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
                       MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
	{
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );

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
                                s += gradientLocalDisplacement[ jcoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) *
                                    gradientLocalDisplacement[ icoor ][ k ][ ig ] * fe.phiDer(i, k, ig ) *
                                    fe.weightDet( ig );
                        }
                    }
                    mat( i, j ) += coef  * s;
                }
            }
	    }
	}
}

// coef * ( \grad \delta u [\grad u^k]^T : \grad v )
void stiff_dergrad_gradbis_2( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
                              MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;

    //
    // blocks (icoor,jcoor) of elmat
    //
    MatrixElemental::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

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
                        s += gradientLocalDisplacement[ l ][ k ][ ig ] * fe.phiDer( i, k, ig ) *
                            fe.phiDer( j, l, ig ) * fe.weightDet( ig );
                }
            }
            mat_tmp( i, j ) = coef * s;
	    }
	}

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
	{
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
	}
}

void stiff_dergrad_gradbis_Tr_2( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
                                 MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;
    //
    // blocks (icoor,jcoor) of elmat
    //
    MatrixElemental::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

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
                        s += gradientLocalDisplacement[ k ][ l ][ ig ] * fe.phiDer( i, k, ig ) *
                            fe.phiDer( j, l, ig ) * fe.weightDet( ig );
                    }
                }
            }
            mat_tmp( i, j ) = coef * s;
	    }
	}

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
	{
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
	}
}



// coef * (  \grad u^k [\grad u]^T \grad u^k : \grad v  )
void stiff_gradgradTr_gradbis_2( Real coef, const boost::multi_array<Real, 3>&  gradientLocalDisplacement,
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

            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor ); // it extracts the (icoor, jcoor) block

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
                                    gradientLocalDisplacement[ jcoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) *
                                    fe.phiDer( j, l, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}

//  coef * (  \grad u [\grad u^k]^T \grad u^k : \grad v  )
void stiff_gradgradTr_gradbis_3( Real coef, const VectorElemental& uk_loc,
                                 MatrixElemental& elmat, const CurrentFE& fe )
{

    Real s;
    // \grad u^k  [\grad u^k]^T  at each quadrature poInt
    //Real guk_gukT[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];
    boost::multi_array<Real, 3> guk_gukT(boost::extents[fe.nbCoor()][fe.nbCoor()][fe.nbQuadPt()]);

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
                            s  += fe.phiDer( i, n, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ] *
                                fe.phiDer( j, n, ig ) * uk_loc.vec() [ j + jcoor * fe.nbFEDof() ] ;
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
    MatrixElemental::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

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
                        s += guk_gukT[ k ][ l ][ ig ] * fe.phiDer( i, k, ig ) *  fe.phiDer( j, l, ig ) * fe.weightDet( ig );
                }
            }
            mat_tmp( i, j ) = coef * s;
	    }
	}

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
	{
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
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
void source_Pvol( Real		coef,
                  const boost::multi_array<Real,3 >& CofFk,
                  const std::vector<Real>& 	Jk,
                  VectorElemental&	elvec,
                  const CurrentFE&	fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        // block (icoor) of elvec
        VectorElemental::vector_view vec =  elvec.block( icoor );
        for( UInt i = 0; i < fe.nbFEDof(); ++i )
	    {
            s = 0.0;
            for( UInt k = 0; k < nDimensions; ++k )
	    	{
                for( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += ( Jk[ig]*Jk[ig] - Jk[ig] + log( Jk[ig] ) )*( 1/Jk[ig] )*
                        CofFk[icoor][ k][ ig]*fe.phiDer(i, k, ig)*fe.weightDet(ig);
                }
	    	}
            vec(i) += coef * s;
	    }
	}
}

//! JACOBIAN MATRIX -------------------------------------------------------------------------------
//! Volumetric part--------------------------------------------------------------------------------

//! 1. Jacobian matrix: Int { 1/2 * coef * ( 2 - 1/J + 1/J^2 ) * ( CofF : \nabla \delta ) (CofF : \nabla v) }
void stiff_Jac_Pvol_1term( Real 	 	 coef,
                           const boost::multi_array<Real,3 >& CofFk,
                           const std::vector<Real>& 	 Jk,
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
                            for( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += ( 2.0 - ( 1/Jk[ig] ) + ( 1/( Jk[ig]*Jk[ig] ) ) ) *
                                    CofFk[ jcoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) *
                                    CofFk[ icoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) *
                                    fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += s * coef;
                }
	    	}
	    }
	}
}


//! 2. Stiffness matrix: Int { 1/2 * coef * ( 1/J - 1 - log(J)/J^2 ) * ( CofF [\nabla \delta]^t CofF ) : \nabla v }
void stiff_Jac_Pvol_2term( Real 		  coef,
                           const boost::multi_array<Real,3 >&  CofFk,
                           const std::vector<Real>& 	  Jk,
                           MatrixElemental&  elmat,
                           const CurrentFE&  fe )
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
                                s +=( ( 1/Jk[ig] ) - 1. - ( 1/( Jk[ig]*Jk[ig] ) ) * log( Jk[ig] ) )*
                                    CofFk[ icoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) *
                                    CofFk[ jcoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += s * coef;
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
void source_P1iso_NH( Real 	    coef,
                      const boost::multi_array<Real,3 >& CofFk,
                      const boost::multi_array<Real,3 >& Fk,
                      const std::vector<Real>&   Jk,
                      const std::vector<Real>&   Ic_isok ,
                      VectorElemental& elvec,
                      const CurrentFE& fe )
{
    Real s1, s2;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        VectorElemental::vector_view vec =  elvec.block( icoor );
        for( UInt i = 0; i < fe.nbFEDof(); ++i )
	    {
            s1 = 0.0; s2 = 0.0;
            for( UInt k = 0; k < nDimensions; ++k )
            {
                for( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s1 +=  pow( Jk[ ig ], (-2.0/3.0) ) * Fk[ icoor ][  k ][ ig ] *
                        fe.phiDer( i, k, ig ) * fe.weightDet( ig );

                    s2 +=  1.0/3.0 * ( Ic_isok[ ig ] * ( 1/Jk[ig] ) ) *
                        CofFk[ icoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                }
            }
            vec( i ) += (s1-s2) * coef;
	    }
	}
}
//! -----------------------------------------------------------------------------------------------






//! Jacobian matrix isochoric part ----------------------------------------------------------------

//! 1. Jacobian matrix : Int { -2/3 * coef * J^(-5/3) *( CofF : \nabla \delta ) ( F : \nabla \v ) }
void stiff_Jac_P1iso_NH_1term( Real coef,
                               const boost::multi_array<Real,3 >& CofFk,
                               const boost::multi_array<Real,3 >& Fk,
                               const std::vector<Real>& Jk ,
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
                            for( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += pow( Jk[ig], -5./3. ) *
                                    Fk[ jcoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) *
                                    CofFk[ icoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += coef * s;
                }
	    	}
	    }
	}
}



//! 2. Stiffness matrix: Int { 2/9 * coef * ( Ic_iso / J^2 )( CofF : \nabla \delta ) ( CofF : \nabla \v ) }
void stiff_Jac_P1iso_NH_2term( Real coef,
                               const boost::multi_array<Real,3 >& CofFk,
                               const std::vector<Real>& Jk ,
                               const std::vector<Real>& Ic_isok,
                               MatrixElemental& elmat,
                               const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor <  nDimensions; ++icoor )
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
                                s += ( 1/(Jk[ig]*Jk[ig] ) ) *  Ic_isok[ig] *
                                    CofFk[ jcoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) *
                                    CofFk[ icoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += coef * s;
                }
	    	}
	    }
	}
}



//! 3. Stiffness matrix : Int { coef * J^(-2/3) (\nabla \delta : \nabla \v)}
void stiff_Jac_P1iso_NH_3term( Real 	     coef,
                               const std::vector<Real>&   Jk,
                               MatrixElemental& elmat,
                               const CurrentFE& fe )
{
    Real s;

    //! assembling diagonal block
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
                    s += pow( Jk[ig], -2./3.) * fe.phiDer( i, k, ig ) *
                        fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                }
	    	}
            mat_tmp( i, j ) = coef * s;
	    }
	}

    for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        //! copy of diagonal block
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
	}
}

//! 4. Stiffness matrix : Int { -2/3 * coef * J^(-5/3) ( F : \nabla \delta ) ( CofF : \nabla \v ) }
void stiff_Jac_P1iso_NH_4term( Real coef,
                               const boost::multi_array<Real,3 >& CofFk,
                               const boost::multi_array<Real,3 >& Fk,
                               const std::vector<Real>& Jk ,
                               MatrixElemental& elmat,
                               const CurrentFE& fe )
{
    Real s;

    for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	{
        for( UInt jcoor = 0; jcoor < nDimensions; ++jcoor )
	    {
            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
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
                                s += pow( Jk[ig], -5./3. ) *
                                    Fk[ icoor ][ k ][ ig ]  * fe.phiDer( i, k, ig ) *
                                    CofFk[ jcoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += coef * s;
                }
	    	}
	    }
	}
}



//! 5. Stiffness matrix : Int { 1/3 * coef * J^(-2) * Ic_iso * (CofF [\nabla \delta]^t CofF ) : \nabla \v }
void stiff_Jac_P1iso_NH_5term( Real coef,
                               const boost::multi_array<Real,3 >& CofFk,
                               const std::vector<Real>& Jk ,
                               const std::vector<Real>& Ic_isok,
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
                                s += ( 1/( Jk[ig]*Jk[ig] ) ) * Ic_isok[ig] *
                                    CofFk[ icoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) *
                                    CofFk[ jcoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += s * coef;
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
void  source_P1iso_Exp( Real             coef,
                        Real             coefExp,
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
                    s += exp( coefExp * ( Ic_isok[ ig ] - 3.0 ) ) *
                        (pow( Jk[ ig ], (-2.0/3.0) ) * Fk[ icoor ][  k ][ ig ] -
                         1.0/3.0 * ( 1/Jk[ ig ] ) * Ic_isok[ ig ] *
                         CofFk[ icoor ][ k ][ ig ] )* fe.phiDer( i, k, ig ) * fe.weightDet( ig );

                }
            }
            vec( i ) += s * coef;
	    }
	}
}

//! Jacobian matrix isochoric part ----------------------------------------------------------------

//! 1. Stiffness term : Int { - 2/3 *coef * J^(-5/3) * exp( coefExp*( Ic_iso - 3) )* ( 1. + coefExp * Ic_iso ) * ( CofF : \nabla \delta ) ( F : \nabla \v ) }
void  stiff_Jac_P1iso_Exp_1term( Real             coef,
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
                                s += pow( Jk[ig], -5./3. ) * exp( coefExp*( Ic_isok[ig] - 3 ) ) *
                                    ( 1. + coefExp * Ic_isok[ig] ) *
                                    CofFk[ icoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) *
                                    Fk[ jcoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) +=  coef * s;
                }
            }
	    }
	}
}

//! 2. Stiffness term : Int { 2 * coef * coefExp * J^(-4/3) * exp( coefExp*( Ic_iso - 3) ) * ( F : \nabla \delta ) ( F : \nabla \v )}
void  stiff_Jac_P1iso_Exp_2term( Real             coef,
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
                                s += pow( Jk[ig], -4.0/3.0 ) * exp( coefExp*(  Ic_isok[ig] -3  ) ) *
                                    Fk[ jcoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) *
                                    Fk[ icoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                            }
                        }
                    }
                    mat( i, j ) += coef * s;
                }
            }
	    }
	}
}

//! 3. Stiffness term: Int { 2.0/9.0 * coef * J^-2 * Ic_iso * exp( coefExp*( Ic_iso - 3) ) * ( 1. + coefExp * Ic_iso )( CofF : \nabla \delta ) ( CofF : \nabla \v )}
void  stiff_Jac_P1iso_Exp_3term( Real coef, Real  coefExp,
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
                                s += ( 1/( Jk[ig]*Jk[ig] ) ) * exp( coefExp*( Ic_isok[ig] - 3 ) ) *
                                    ( 1. + coefExp * Ic_isok[ig] )* Ic_isok[ig] *
                                    CofFk[ jcoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) *
                                    CofFk[ icoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) * fe.weightDet( ig );

                            }
                        }
                    }
                    mat( i, j ) +=  coef * s;
                }
            }
	    }
	}
}

//! 4. Stiffness term: Int { -2.0/3.0 * coef * J^(-5/3) * exp( coefExp*( Ic_iso - 3) ) * ( 1. + coefExp * Ic_iso )( F : \nabla \delta ) ( CofF : \nabla \v ) }
void  stiff_Jac_P1iso_Exp_4term( Real coef, Real  coefExp,
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
                                s += pow( Jk[ig], -5./3. ) * exp( coefExp*( Ic_isok[ig] - 3  ) ) *
                                    ( 1. + coefExp * Ic_isok[ig] ) *
                                    Fk[ icoor ][ k ][ ig ]  * fe.phiDer( i, k, ig ) *
                                    CofFk[ jcoor ][ l ][ ig ] *  fe.phiDer( j, l, ig ) * fe.weightDet( ig );
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
void  stiff_Jac_P1iso_Exp_5term( Real             coef,
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
                    s += pow(Jk[ig], -2.0/3.0) * exp( coefExp*( Ic_isok[ig] -3  ) ) *
                        fe.phiDer( i, k, ig ) *  fe.phiDer( j, k, ig ) * fe.weightDet( ig );
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

//! 6. Stiffness term : Int { 1.0/3.0 * coef * J^(-2) * Ic_iso *  exp(coefExp( Ic_iso - 3)) * (CofF [\nabla \delta]^t CofF ) : \nabla \v }
void  stiff_Jac_P1iso_Exp_6term( Real             coef,
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
                                s += ( 1/( Jk[ig]*Jk[ig] ) ) * Ic_isok[ig] *
                                    exp( coefExp*( Ic_isok[ig] -3  ) ) *
                                    CofFk[ icoor ][ l ][ ig ] * fe.phiDer( i, k, ig ) *
                                    CofFk[ jcoor ][ k ][ ig ] * fe.phiDer( j, l, ig ) * fe.weightDet( ig );
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
//! END OF EXPONENTIAL MODEL
//! ***********************************************************************************************

} //! End namespace AssemblyElementalStructure

} //! End namespace LifeV

#endif
