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
 *  @file
 *  @brief  This class has been created only to split the methods for structural problems fr *          om the others methods of AssemblyElemental. It contains all the methods are priv *          ate.
 *          All the methods have been developed to implement the materials which are in Life *          V at the moment of writing: linear elastic, Venant-Kirchhoff, neohookean and exp *          onential.
 *
 *  @version 1.0
 *  @date 28-08-2011
 *  @author Paolo Tricerri
 *
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _STRUCTURALASSEMBLER_CPP_
#define _STRUCTURALASSEMBLER_CPP_

#include<life/lifesolver/StructuralAssembler.hpp>

namespace LifeV
{

//==================================
// Constructor & Decostructor
//==================================
StructuralAssembler::StructuralAssembler()
{}

StructuralAssembler::~StructuralAssembler()
{}

//==================================
// Public Methods
//==================================
  
//! Methods for the linear elastic model

void StructuralAssembler::stiff_strain( Real coef, MatrixElemental& elmat, const CurrentFE& fe )
{
    double s;
    double tmp = coef * 0.5;

    MatrixElemental::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0;
            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
                    s += fe.phiDer( i, icoor, ig ) * fe.phiDer( j, icoor, ig ) * fe.weightDet( ig );
            mat_tmp( i, j ) = tmp * s;
        }
    }
    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
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
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                        s += fe.phiDer( i, jcoor, ig ) * fe.phiDer( j, icoor, ig ) * fe.weightDet( ig );
                    mat( i, j ) += tmp * s;
                }
            }
        }
    }
}

  
void StructuralAssembler::stiff_div( Real coef, MatrixElemental& elmat, const CurrentFE& fe )
{
    double s;

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
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                        s += fe.phiDer( i, icoor, ig ) * fe.phiDer( j, jcoor, ig ) * fe.weightDet( ig );
                    mat( i, j ) += coef * s;
                }
            }
        }
    }
}
// End of methods for linear elastic model

//! Methods for St. Venant Kirchhoff model
//! Methods for the stiffness matrix 
void StructuralAssembler::stiff_derdiv( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    Real guk[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k at each quadrature point
    Real s;

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
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
                            s += fe.phiDer( i, icoor, ig ) * guk[ jcoor ][ k ][ ig ] * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                    mat( i, j ) += coef * s;
                }
            }
        }
    }
}

void StructuralAssembler::stiff_dergradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real guk[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k at each quadrature point


    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
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
                            s += fe.phiDer( i, k, ig ) * guk[ jcoor ][ icoor ][ ig ] * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                    mat( i, j ) += coef * s;
                }
            }
        }
    }
}

void StructuralAssembler::stiff_divgrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real duk[ fe.nbQuadPt() ];      // definizione del vettore che conterra \div u^k at each quadrature point


    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        s=0;
        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {
            // s = 0.0; Alessandro
            for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                s += fe.phiDer( i, icoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; // costruzione di \div u^k at a quadrature point

            // duk[ ig ] = s; Alessandro
        }// chiude il ciclo su icoor
        duk[ ig ] = s;
    }// chiude il ciclo su ig

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

void StructuralAssembler::stiff_gradgrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s,s1;
    Real gguk[ fe.nbQuadPt() ]; //    (\grad u_k : \grad u_k) at each quadrature point

    // loop on quadrature points
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
        }// chiude il ciclo su icoor
        gguk[ ig ] = s;
    }// chiude il ciclo su ig

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

void StructuralAssembler::stiff_dergrad_gradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real guk[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k at each quadrature point


    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
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
                            s += guk[ icoor ][ jcoor ][ ig ] * fe.phiDer( i, k, ig ) *  fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                    }
                    mat( i, j ) += coef * s;
                }
            }
        }
    }
}

void StructuralAssembler::stiff_dergrad_gradbis_Tr( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real guk[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k at each quadrature point


    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
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
                            s += guk[ icoor ][ k ][ ig ]  * fe.phiDer( j, k, ig ) * fe.phiDer( i, jcoor, ig ) * fe.weightDet( ig );
                    }
                    mat( i, j ) += coef * s;
                }
            }
        }
    }
}

void StructuralAssembler::stiff_gradgradTr_gradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real guk_gukT[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k  [\grad u^k]^T  at each quadrature point

    // loop on quadrature points
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
                            s  += fe.phiDer( i, n, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ] * fe.phiDer( j, n, ig ) * uk_loc.vec() [ j + jcoor * fe.nbFEDof() ] ; // \grad u^k  [\grad u^k]^T  at each quadrature point
                    }
                }
                guk_gukT[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor ); // estrae il blocco (icoor, jcoor)

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            s += fe.phiDer( i, k, ig ) * guk_gukT[ icoor ][ jcoor ][ ig ] * fe.phiDer( j, k, ig ) * fe.weightDet( ig );
                    mat( i, j ) += coef * s;
                }
            }
        }
    }
}

// End of methods for the stiffness matrix (St. Venant-Kirchhoff material)
//! Methods for the jacobian

void StructuralAssembler::stiff_dergrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real guk[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k at each quadrature point


    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
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
                        {
                            s += fe.phiDer( i, k, ig ) * ( guk[ jcoor ][ k ][ ig ] * fe.phiDer( j, icoor, ig )
                                                           + guk[ jcoor ][ icoor ][ ig ] * fe.phiDer( j, k, ig ) ) * fe.weightDet( ig );
                        }
                    mat( i, j ) += coef * s;
                }
            }
        }
    }
}

void StructuralAssembler::stiff_divgrad_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real guk[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k at each quadrature point

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {
            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }

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
                            s += fe.phiDer( j, jcoor, ig ) * guk[ icoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
                    mat( i, j ) += coef * s;
                }
            }
        }
    }
}

void StructuralAssembler::stiff_gradgrad_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real guk[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k at each quadrature point

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {
            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
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
                    s = 0.0;
                    for ( UInt k = 0; k < fe.nbCoor(); ++k )
                    {
                        for ( UInt l = 0; l < fe.nbCoor(); ++l )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                                s += guk[ jcoor ][ l ][ ig ] * fe.phiDer( j, l, ig ) * guk[ icoor ][ k ][ ig ] * fe.phiDer(i, k, ig ) * fe.weightDet( ig ); // il ciclo sui nodi di quadratura
                        }       																																							 // fa l'integrale
                    }
                    mat( i, j ) += coef  * s;
                }
            }
        }
    }
}

void StructuralAssembler::stiff_dergrad_gradbis_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real guk[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k at each quadrature point


    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
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
                        s += guk[ l ][ k ][ ig ] * fe.phiDer( i, k, ig ) *  fe.phiDer( j, l, ig ) * fe.weightDet( ig );
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

void StructuralAssembler::stiff_dergrad_gradbis_Tr_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real guk[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k at each quadrature point


    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
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
                        s += guk[ k ][ l ][ ig ] * fe.phiDer( i, k, ig ) *  fe.phiDer( j, l, ig ) * fe.weightDet( ig );
                }
            }
            mat_tmp( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor ) // copia del blocco sulla diagonale
    {
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
    }
}

void StructuralAssembler::stiff_gradgradTr_gradbis_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real guk[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k at each quadrature point


    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }

    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor ); // estrae il blocco (icoor, jcoor)

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
                                s += guk[ icoor ][ l ][ ig ] *guk[ jcoor ][ k ][ ig ] * fe.phiDer( i, k, ig ) *  fe.phiDer( j, l, ig ) * fe.weightDet( ig );
                        }
                    }
                    mat( i, j ) += coef * s;
                }
            }
        }
    }
}

void StructuralAssembler::stiff_gradgradTr_gradbis_3( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    Real guk_gukT[ fe.nbCoor() ][ fe.nbCoor() ][ fe.nbQuadPt() ];      // \grad u^k  [\grad u^k]^T  at each quadrature point
    // attenzione in questa funzione si deve usare il trasposto
    // loop on quadrature points                                                // (\grad u^k  [\grad u^k]^T )^T
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
                            s  += fe.phiDer( i, n, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ] * fe.phiDer( j, n, ig ) * uk_loc.vec() [ j + jcoor * fe.nbFEDof() ] ; // \grad u^k  [\grad u^k]^T  at each quadrature point
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

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor ) // copia del blocco sulla diagonale
    {
        MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
    }
}
// End of St. Venant Kirchhoff model

//! Methods for NeoHookean model
void StructuralAssembler::source_Pvol( Real			coef, 
				       const	KNMK<Real> 	CofFk, 
				       const	KN<Real> 	Jk, 
				       VectorElemental&		elvec, 
				       const CurrentFE&	fe )
{
  double s;
  // ASSERT_PRE( fe.hasFirstDeriv(),
  //	"source_stress needs at least the velocity shape functions first derivatives" );

  for ( int icoor = 0; icoor < nDimensions; ++icoor )
    {
      // block (icoor) of elvec
      VectorElemental::vector_view vec =  elvec.block( icoor ); 
      for ( int i = 0; i < fe.nbFEDof(); ++i )
	{
	  s = 0.0;
	  for ( int k = 0; k < nDimensions; ++k )
	    {
	      for ( int ig = 0; ig < fe.nbQuadPt(); ++ig )
		{
		  s += ( pow( Jk(ig),2 ) - Jk(ig) + log( Jk(ig) ) )*pow( Jk(ig),-1)*
		    CofFk(icoor, k, ig)*fe.phiDer(i, k, ig)*fe.weightDet(ig);  
		  //s +=( Jk( ig ) - 1.  + pow( Jk( ig ), -1 )  *log( Jk( ig ) )  ) * 
		  //CofFk( icoor, k, ig ) * fe.phiDer( i, k, ig )  * fe.weightDet( ig );
		}
	    }
	  vec(i) += coef/2.0 * s;
	}
    }
}

void StructuralAssembler::source_P1iso_NH(Real coef, 
					  const KNMK<Real> CofFk, 
					  const KNMK<Real> Fk,
					  const KN<Real> Jk, 
					  const KN<Real>Ic_isok , 
					  VectorElemental& elvec, 
					  const CurrentFE& fe)
{
    
  //ASSERT_PRE( fe.hasFirstDeriv(),
  //	"source_stress needs at least the velocity shape functions first derivatives" );

  double s1, s2;
  
  for ( int icoor = 0; icoor < nDimensions; ++icoor )
    {
      
      // block (icoor) of elvec
      VectorElemental::vector_view vec =  elvec.block( icoor ); 
      for ( int i = 0; i < fe.nbFEDof(); ++i ){
      
	s1 = 0.0; s2 = 0.0;
	for ( int k = 0; k < nDimensions; ++k )
	  {	 
	    for ( int ig = 0;ig < fe.nbQuadPt(); ++ig )
	      {
	    
		s1 +=  pow( Jk( ig ), (-2.0/3.0) ) * Fk( icoor,  k, ig ) * 
		  fe.phiDer( i, k, ig ) * fe.weightDet( ig );
	    
		s2 +=  1.0/3.0 * ( Ic_isok( ig ) * pow( Jk( ig ), -1.0 ) ) * 
		  CofFk( icoor, k, ig ) * fe.phiDer( i, k, ig ) * fe.weightDet( ig )  ;
	    
	      }
	  }
	vec( i ) += (s1-s2) * coef  ;
      }
    }  
}

void StructuralAssembler::stiff_Jac_Pvol_1term( Real 	 	 coef, 
						const KNMK<Real> CofFk, 
						const KN<Real> 	 Jk, 
						MatrixElemental& elmat, 
						const CurrentFE& fe )
{
  //ASSERT_PRE( fe.hasFirstDeriv(),
  //        "Stiffness dergradbis matrix needs at least the first derivatives" );
    
  double s;
    
  //! blocks (icoor,jcoor) of elmat
  for ( int icoor = 0; icoor < nDimensions; ++icoor )
    {     
      for ( int jcoor = 0; jcoor < nDimensions; ++jcoor )
	{
	  MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor ); 
	  for ( int i = 0; i < fe.nbFEDof(); ++i )
	    {
	      for ( int j = 0; j < fe.nbFEDof(); ++j )
		{
		  s = 0.0;
		  for ( int l = 0; l < nDimensions; ++l )
		    {
		      for ( int k = 0; k < nDimensions; ++k )
			{		
			  for ( int ig = 0;ig < fe.nbQuadPt(); ++ig )
			    {
			      s += ( 2.0 -   pow(Jk(ig), -1.) + pow(Jk(ig), -2.)  ) *
				CofFk( jcoor , l , ig ) * fe.phiDer( j, l, ig ) *  
				CofFk( icoor , k , ig ) * fe.phiDer( i, k, ig ) * 
				fe.weightDet( ig );
			    }
			}
		    }
		  mat( i, j ) += s * coef/2.0 ;
		}
	    }
	}
    }
}

void StructuralAssembler::stiff_Jac_Pvol_2term( Real 		  coef, 
						const KNMK<Real>  CofFk, 
						const KN<Real> 	  Jk,  
						MatrixElemental&  elmat, 
						const CurrentFE&  fe )
{
  // ASSERT_PRE( fe.hasFirstDeriv(),
  //        "Stiffness dergradbis matrix needs at least the first derivatives" );
    
  double s;
    
  // blocks (icoor,jcoor) of elmat
  for ( int icoor = 0; icoor < nDimensions; ++icoor )
    {      
      for ( int jcoor = 0; jcoor < nDimensions; ++jcoor )
	{	
	  MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor ); 
	  for ( int i = 0; i < fe.nbFEDof(); ++i )
	    {	  
	      for ( int j = 0; j < fe.nbFEDof(); ++j )
		{	    
		  s = 0.0;
		  for ( int l = 0; l < nDimensions; ++l )
		    {	      
		      for ( int k = 0; k < nDimensions; ++k )
			{
			  for ( int ig = 0;ig < fe.nbQuadPt(); ++ig )
			    {
			      s +=(  pow(Jk(ig), -1.)   - 1.  - pow(Jk(ig), -2) * log( Jk(ig) )   )*
				CofFk( icoor , l , ig ) * fe.phiDer( j, l, ig ) * 
				CofFk( jcoor , k , ig ) * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
			    }
			}
		    }
		  mat( i, j ) += s * coef /2.0 ;
		}
	    }
	}
    }
}

void StructuralAssembler::stiff_Jac_P1iso_NH_1term( Real coef, 
						    const KNMK<Real> CofFk, 
						    const KNMK<Real> Fk,
						    const KN<Real> Jk , 
						    MatrixElemental& elmat, 
						    const CurrentFE& fe )
{
  // ASSERT_PRE( fe.hasFirstDeriv(),
  //           "Stiffness dergradbis matrix needs at least the first derivatives" );
  
  double s;
    
  // blocks (icoor,jcoor) of elmat
  //
  for ( int icoor = 0; icoor < nDimensions; ++icoor )
    {
      for ( int jcoor = 0; jcoor < nDimensions; ++jcoor )
	{
	  MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor ); 
	
	  for ( int i = 0; i < fe.nbFEDof(); ++i )
	    {
	      for ( int j = 0; j < fe.nbFEDof(); ++j )
		{
		  s = 0.0;	    
		  for ( int l = 0; l < nDimensions; ++l )
		    {
		      for ( int k = 0; k < nDimensions; ++k )
			{		
			  for ( int ig = 0;ig < fe.nbQuadPt(); ++ig )
			    s += pow( Jk(ig), -5./3. ) * 
			      Fk( jcoor , l , ig ) * fe.phiDer( j, l, ig ) * 
			      CofFk( icoor , k , ig ) * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
		
			}
		    }
		  mat( i, j ) += -2.0/3.0 * coef * s;
		}
	    }
	}
    }
}

void StructuralAssembler::stiff_Jac_P1iso_NH_2term( Real coef, 
						    const KNMK<Real> CofFk,
						    const KN<Real> Jk , 
						    const KN<Real> Ic_isok, 
						    MatrixElemental& elmat, 
						    const CurrentFE& fe )
{
  //ASSERT_PRE( fe.hasFirstDeriv(),
  //          "Stiffness dergradbis matrix needs at least the first derivatives" );
    
  double s;
  
  // blocks (icoor,jcoor) of elmat
  for ( int icoor = 0; icoor < nDimensions; ++icoor )
    {
      for ( int jcoor = 0; jcoor < nDimensions; ++jcoor )
	{
	  MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor ); 
	  for ( int i = 0; i < fe.nbFEDof(); ++i )
	    {
	      for ( int j = 0; j < fe.nbFEDof(); ++j )
		{
		  s = 0.0;
		  for ( int l = 0; l < nDimensions; ++l )
		    {
		      for ( int k = 0; k < nDimensions; ++k )
			{
			  for ( int ig = 0;ig < fe.nbQuadPt(); ++ig )
			    {
			      s += pow( Jk(ig), -2. ) *  Ic_isok(ig)  *
				CofFk( jcoor , l , ig )  * fe.phiDer( j, l, ig ) *
				CofFk( icoor , k , ig ) *  fe.phiDer( i, k, ig ) * fe.weightDet( ig );
			    }
			}
		    }
		  mat( i, j ) += 2./9. * coef * s;
		}
	    }
	}
    }
}

void StructuralAssembler::stiff_Jac_P1iso_NH_3term( Real coef, 
						    const KN<Real> Jk , 
						    MatrixElemental& elmat, 
						    const CurrentFE& fe )
{
  // ASSERT_PRE( fe.hasFirstDeriv(),
  //	"Stiffness dergradbis matrix needs at least the first derivatives" );
    
  double s;
    
  // assembling diagonal block
  MatrixElemental::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    
  for ( int i = 0; i < fe.nbFEDof(); ++i )
    {
      for ( int j = 0; j < fe.nbFEDof(); ++j )
	{
	  s = 0.0;
	  for ( int k = 0; k < nDimensions; ++k )
	    {
	      for ( int ig = 0;ig < fe.nbQuadPt(); ++ig )
		s += pow( Jk(ig), -2./3.) * fe.phiDer( i, k, ig ) *  
		  fe.phiDer( j, k, ig ) * fe.weightDet( ig );
	    }
	  mat_tmp( i, j ) = coef * s;
	}
    }
    
  for ( int icoor = 0; icoor < nDimensions; ++icoor ) 
    {// copy of diagonal block
      MatrixElemental::matrix_view mat = elmat.block( icoor, icoor );
      mat += mat_tmp;
    }  
}

void StructuralAssembler::stiff_Jac_P1iso_NH_4term( Real coef, 
						    const KNMK<Real> CofFk, 
						    const KNMK<Real> Fk,
						    const KN<Real> Jk , 
						    MatrixElemental& elmat, 
						    const CurrentFE& fe )
{
  //ASSERT_PRE( fe.hasFirstDeriv(),
  //         "Stiffness dergradbis matrix needs at least the first derivatives" );
    
  double s;
    
  // blocks (icoor,jcoor) of elmat
  for ( int icoor = 0; icoor < nDimensions; ++icoor )
    {
      for ( int jcoor = 0; jcoor < nDimensions; ++jcoor )
	{
	  MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor ); 	
	  for ( int i = 0; i < fe.nbFEDof(); ++i )
	    {	  
	      for ( int j = 0; j < fe.nbFEDof(); ++j )
		{
		  s = 0.0;
		  for ( int l = 0; l < nDimensions; ++l )
		    {	      
		      for ( int k = 0; k < nDimensions; ++k )
			{		
			  for ( int ig = 0;ig < fe.nbQuadPt(); ++ig )
			    s += pow( Jk(ig), -5./3. ) * 
			      Fk( icoor , k , ig )  * fe.phiDer( i, k, ig ) * 
			      CofFk( jcoor , l , ig ) *  fe.phiDer( j, l, ig ) * fe.weightDet( ig );
			}
		    }
		  mat( i, j ) += - 2./3. * coef * s;
		}
	    }
	}
    }
}

void StructuralAssembler::stiff_Jac_P1iso_NH_5term( Real coef,  
						    const KNMK<Real> CofFk, 
						    const KN<Real> Jk ,
						    const KN<Real> Ic_isok, 
						    MatrixElemental& elmat, 
						    const CurrentFE& fe )
{
  //ASSERT_PRE( fe.hasFirstDeriv(),
  //      "Stiffness dergradbis matrix needs at least the first derivatives" );
    
  double s;
    
  // blocks (icoor,jcoor) of elmat
  for ( int icoor = 0; icoor < nDimensions; ++icoor )
    {
      for ( int jcoor = 0; jcoor < nDimensions; ++jcoor )
	{
	  MatrixElemental::matrix_view mat = elmat.block( icoor, jcoor ); 
	
	  for ( int i = 0; i < fe.nbFEDof(); ++i )
	    {
	      for ( int j = 0; j < fe.nbFEDof(); ++j )
		{
		  s = 0.0;
		  for ( int l = 0; l < nDimensions; ++l )
		    {
		      for ( int k = 0; k < nDimensions; ++k )
			{
			  for ( int ig = 0;ig < fe.nbQuadPt(); ++ig )
			    {
			      s += pow(Jk(ig), -2.) *  Ic_isok(ig) *
				CofFk( icoor , l , ig ) * fe.phiDer( j, l, ig ) * 
				CofFk( jcoor , k , ig ) * fe.phiDer( i, k, ig ) * fe.weightDet( ig );
			    }
			}
		    }
		  mat( i, j ) +=  1./3. * s * coef ;
		}
	    }
	}
    }
}
// End of Neohookean model

}
#endif /* _STRUCTURALASSEMBLER_CPP_ */
