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
    @brief File containing the procedures for the local assembly of the differential operators

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef ELEMOPER_CPP
#define ELEMOPER_CPP 1

#include <life/lifefem/AssemblyElemental.hpp>

namespace LifeV
{

namespace ElemOper
{
void mass(ElemMat& localMass,
          const CurrentFE& massCFE,
          const Real& coefficient,
          const UInt& fieldDim)
{
    const UInt nbFEDof(massCFE.nbFEDof());
    const UInt nbQuadPt(massCFE.nbQuadPt());
    Real localValue(0);

    // Assemble the local mass
    for (UInt iterFDim(0); iterFDim<fieldDim; ++iterFDim)
    {
        // Extract the view of the matrix
        ElemMat::matrix_view localView = localMass.block(iterFDim,iterFDim);

        // Loop over the basis functions
        for (UInt iDof(0); iDof < nbFEDof ; ++iDof)
        {
            // Build the local matrix only where needed:
            // Lower triangular + diagonal parts
            for (UInt jDof(0); jDof <= iDof; ++jDof)
            {
                localValue = 0.0;

                //Loop on the quadrature nodes
                for (UInt iQuadPt(0); iQuadPt < nbQuadPt; ++iQuadPt)
                {
                    localValue += massCFE.phi(iDof,iQuadPt)
                                  * massCFE.phi(jDof,iQuadPt)
                                  * massCFE.wDetJacobian(iQuadPt);
                }

                localValue*=coefficient;

                // Add on the local matrix
                localView(iDof,jDof)+=localValue;

                if (iDof!=jDof)
                {
                    localView(jDof,iDof)+=localValue;
                }
            }
        }
    }
}


void stiffness(ElemMat& localStiff,
               const CurrentFE& stiffCFE,
               const Real& coefficient,
               const UInt& fieldDim)
{
    const UInt nbFEDof(stiffCFE.nbFEDof());
    const UInt nbQuadPt(stiffCFE.nbQuadPt());
    Real localValue(0);

    // Assemble the local diffusion
    for (UInt iterFDim(0); iterFDim<fieldDim; ++iterFDim)
    {
        // Extract the view of the matrix
        ElemMat::matrix_view localView = localStiff.block(iterFDim,iterFDim);

        // Loop over the basis functions
        for (UInt iDof(0); iDof < nbFEDof ; ++iDof)
        {
            // Build the local matrix only where needed:
            // Lower triangular + diagonal parts
            for (UInt jDof(0); jDof <= iDof; ++jDof)
            {
                localValue = 0.0;

                //Loop on the quadrature nodes
                for (UInt iQuadPt(0); iQuadPt < nbQuadPt; ++iQuadPt)
                {
                    for (UInt iDim(0); iDim<3; ++iDim)
                    {
                        localValue += stiffCFE.dphi(iDof,iDim,iQuadPt)
                                      * stiffCFE.dphi(jDof,iDim,iQuadPt)
                                      * stiffCFE.wDetJacobian(iQuadPt);
                    }
                }

                localValue*=coefficient;

                // Add on the local matrix
                localView(iDof,jDof)+=localValue;

                if (iDof != jDof)
                {
                    localView(jDof,iDof)+=localValue;
                }
            }
        }
    }
}
}

//
//----------------------------------------------------------------------
//                      Element matrix operator
//----------------------------------------------------------------------
//Real coef(t,x,y,z,u)
/*
  Mass matrix: \int coef(t,x,y,z,u) v_i v_j
*/

//
// coeff*Mass
//
void mass( Real coef, ElemMat& elmat, const CurrentFE& fe,
           int iblock, int jblock )
/*
  Mass matrix: \int v_i v_j
*/
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt i, ig;
    int iloc, jloc;
    Real s, coef_s;
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += fe.phi( iloc, ig ) * fe.phi( iloc, ig ) * fe.weightDet( ig );
        }
        mat( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
            s += fe.phi( iloc, ig ) * fe.phi( jloc, ig ) * fe.weightDet( ig );
        coef_s = coef * s;
        mat( iloc, jloc ) += coef_s;
        mat( jloc, iloc ) += coef_s;
    }
}
//
// coeff*Mass
//
void mass( Real coef, ElemMat& elmat, const CurrentFE& fe,
           int iblock, int jblock, UInt nb )
/*
  Mass matrix: \int v_i v_j (nb blocks on the diagonal, nb>1)
*/
{
    Matrix mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    UInt i, ig;
    int iloc, jloc;
    Real s, coef_s;
    mat_tmp = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += fe.phi( iloc, ig ) * fe.phi( iloc, ig ) * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
            s += fe.phi( iloc, ig ) * fe.phi( jloc, ig ) * fe.weightDet( ig );
        coef_s = coef * s;
        mat_tmp( iloc, jloc ) += coef_s;
        mat_tmp( jloc, iloc ) += coef_s;
    }
    // copy on the components
    for ( UInt icomp = 0; icomp < nb; icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst( i );
            jloc = fe.patternSecond( i );
            mat_icomp( iloc, jloc ) += mat_tmp( iloc, jloc );
            mat_icomp( jloc, iloc ) += mat_tmp( jloc, iloc );
        }
    }
}

//
// coeff[q]*Mass
//
void mass( const std::vector<Real>& coef, ElemMat& elmat, const CurrentFE& fe,
           int iblock, int jblock, UInt nb )
/*
  Mass matrix: \int v_i v_j (nb blocks on the diagonal, nb>1)
*/
{
    Matrix mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    UInt i, ig;
    int iloc, jloc;
    Real s;//, coef_s;
    mat_tmp = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += coef[ig] * fe.phi( iloc, ig ) * fe.phi( iloc, ig ) * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) = s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
            s += coef[ig]*fe.phi( iloc, ig ) * fe.phi( jloc, ig ) * fe.weightDet( ig );
        mat_tmp( iloc, jloc ) += s;
        mat_tmp( jloc, iloc ) += s;
    }
    // copy on the components
    for ( UInt icomp = 0; icomp < nb; icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst( i );
            jloc = fe.patternSecond( i );
            mat_icomp( iloc, jloc ) += mat_tmp( iloc, jloc );
            mat_icomp( jloc, iloc ) += mat_tmp( jloc, iloc );
        }
    }
}






void stiff_divgrad( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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

    ElemMat::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

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
        ElemMat::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
    }

}

// Stiffness matrix: coef * ( (\div u) \grad u_k : \grad v  ) controllato!!!
void stiff_divgrad_2( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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
            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

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

// Stiffness matrix: coef * ( \grad u_k : \grad u_k) *( \grad u : \grad v  ) controllato!!!
void stiff_gradgrad( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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

    ElemMat::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

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
        ElemMat::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
    }
}

// Stiffness matrix: coef * ( \grad u_k : \grad u) *( \grad u_k : \grad v  ) controllato!!!
void stiff_gradgrad_2( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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
            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

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

// Stiffness matrix: coef * ( \grad d^k \grad d : \grad v  )controllato!!!
void stiff_dergrad_gradbis( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

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

// Stiffness matrix: coef * ( \grad u^k [\grad d]^T : \grad v  ) controllato!!!
void stiff_dergrad_gradbis_Tr( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

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

// Stiffness matrix: coef * ( \grad \delta d \grad d^k : \grad v  ) controllato!!!
void stiff_dergrad_gradbis_2( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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
    ElemMat::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

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
        ElemMat::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
    }
}

// Stiffness matrix: coef * ( \grad \delta d [\grad d^k]^T : \grad v  )
void stiff_dergrad_gradbis_Tr_2( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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
    ElemMat::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

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
        ElemMat::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
    }
}

//  Stiffness matrix: coef * (  \grad u^k [\grad u^k]^T \grad u : \grad v  ) controllato!!!
void stiff_gradgradTr_gradbis( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor ); // estrae il blocco (icoor, jcoor)

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

//  Stiffness matrix: coef * (  \grad u^k [\grad u]^T \grad u^k : \grad v  )controllato!!!
void stiff_gradgradTr_gradbis_2( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor ); // estrae il blocco (icoor, jcoor)

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

//  Stiffness matrix: coef * (  \grad u [\grad u^k]^T \grad u^k : \grad v  )controllato!!!
void stiff_gradgradTr_gradbis_3( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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
    ElemMat::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

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
        ElemMat::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
    }
}







void ipstab_grad( const Real         coef,
                  ElemMat&           elmat,
                  const CurrentFE&   fe1,
                  const CurrentFE&   fe2,
                  const CurrentBdFE& bdfe,
                  int iblock, int jblock )
{
    /*
      Interior penalty stabilization: coef*\int_{face} grad u1_i . grad v1_j
    */

    ElemMat::matrix_view mat = elmat.block( iblock, jblock );



    Real sum, sum1, sum2;
    UInt icoor,jcoor,i,j;
    UInt ig;
    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];
    Real phid1[ fe1.nbFEDof() ][ fe1.nbCoor() ][ bdfe.nbQuadPt() ];
    Real phid2[ fe2.nbFEDof() ][ fe2.nbCoor() ][ bdfe.nbQuadPt() ];
    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac( icoor, jcoor, 0 ) * drp2[ jcoor ];
                }
                phid1[ i ][ icoor ][ ig ] = sum1;
                phid2[ i ][ icoor ][ ig ] = sum2;
            }
        }
    }



    // Loop on rows
    for ( i = 0; i < fe1.nbFEDof(); ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbFEDof(); ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt() ; ++ig )
                    sum += phid1[ i ][ icoor ][ ig ] * phid2[ j ][ icoor ][ ig ] * bdfe.weightMeas( ig );
            mat( i, j ) += coef * sum;
        }
    }

}






void ipstab_grad( const Real         coef,
                  ElemMat&           elmat,
                  const CurrentFE&   fe1,
                  const CurrentFE&   fe2,
                  const CurrentBdFE& bdfe,
                  int iblock, int jblock,
                  int nb )
{
    /*
      Interior penalty stabilization: coef*\int_{face} grad u1_i . grad v1_j
    */


    ElemMat::matrix_type mat_tmp( fe1.nbFEDof(), fe2.nbFEDof() );


    Real sum, sum1, sum2;
    UInt icoor,jcoor,i,j;
    UInt ig;
    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];
    Real phid1[ fe1.nbFEDof() ][ fe1.nbCoor() ][ bdfe.nbQuadPt() ];
    Real phid2[ fe2.nbFEDof() ][ fe2.nbCoor() ][ bdfe.nbQuadPt() ];
    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac( icoor, jcoor, 0 ) * drp2[ jcoor ];
                }
                phid1[ i ][ icoor ][ ig ] = sum1;
                phid2[ i ][ icoor ][ ig ] = sum2;
            }
        }
    }


    // Loop on rows
    for ( i = 0; i < fe1.nbFEDof(); ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbFEDof(); ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt() ; ++ig )
                    sum += phid1[ i ][ icoor ][ ig ] * phid2[ j ][ icoor ][ ig ] * bdfe.weightMeas( ig );
            mat_tmp( i, j ) = coef * sum;
        }
    }


    // copy on the components
    for ( int icomp = 0; icomp < nb; icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        mat_icomp += mat_tmp;
    }
}





void ipstab_bgrad( const Real         coef,
                   ElemMat&           elmat,
                   const CurrentFE&   fe1,
                   const CurrentFE&   fe2,
                   const ElemVec&     beta,
                   const CurrentBdFE& bdfe,
                   int iblock, int jblock,
                   int nb )
{
    /*
      Interior penalty stabilization: coef*\int_{face} (\beta1 . grad u1_i) . (\beta2 . grad v2_j)
    */

    ElemMat::matrix_type mat_tmp( fe1.nbFEDof(), fe2.nbFEDof() );

    Real sum, sum1, sum2;
    UInt i,j;
    UInt icoor,jcoor;
    UInt ig;

    //
    // convection velocity \beta on the boundary quadrature points
    //
    Real b[ fe1.nbCoor() ][ bdfe.nbQuadPt() ];

    for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
    {
        for ( ig = 0; ig < bdfe.nbQuadPt(); ig++ )
        {
            sum = 0;
            for ( i = 0; i < bdfe.nbNode(); ++i )
            {
                sum += bdfe.phi( i, ig ) * beta.vec() [ icoor * bdfe.nbCoor() + i ];
            }
            b[ icoor ][ ig ] = sum;
        }
    }


    //
    // shape function first derivaties on the boundary quadrature points
    //
    // this should be improved!!!
    //

    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];
    Real phid1[ fe1.nbFEDof() ][ fe1.nbCoor() ][ bdfe.nbQuadPt() ];
    Real phid2[ fe2.nbFEDof() ][ fe2.nbCoor() ][ bdfe.nbQuadPt() ];
    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac( icoor, jcoor, 0 ) * drp2[ jcoor ];
                }
                phid1[ i ][ icoor ][ ig ] = sum1;
                phid2[ i ][ icoor ][ ig ] = sum2;
            }
        }
    }

    // Loop on rows
    for ( i = 0; i < fe1.nbFEDof(); ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbFEDof(); ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
                for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
                    for ( ig = 0; ig < bdfe.nbQuadPt(); ig++ )
                        sum += phid1[ i ][ icoor ][ ig ]*phid2[ j ][ jcoor ][ ig ]
                               *b[ icoor ][ ig ]*b[ jcoor ][ ig ]
                               *bdfe.weightMeas( ig );
            mat_tmp( i, j ) = coef * sum;
        }
    }

    // copy on the components
    for ( int icomp = 0; icomp < nb; icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        mat_icomp += mat_tmp;
    }

}





void ipstab_div( const Real coef, ElemMat& elmat, const CurrentFE& fe1, const CurrentFE& fe2,
                 const CurrentBdFE& bdfe, int iblock, int jblock )
{
    /*
      Interior penalty stabilization: coef*\int_{face} div u . div v
    */

    Real sum, sum1, sum2;
    UInt i,j,icoor,jcoor;
    UInt ig;
    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];
    Real phid1[ fe1.nbFEDof() ][ fe1.nbCoor() ][ bdfe.nbQuadPt() ];
    Real phid2[ fe2.nbFEDof() ][ fe2.nbCoor() ][ bdfe.nbQuadPt() ];
    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac( icoor, jcoor, 0 ) * drp2[ jcoor ];
                }
                phid1[ i ][ icoor ][ ig ] = sum1;
                phid2[ i ][ icoor ][ ig ] = sum2;
            }
        }
    }

    for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
    {
        for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
        {
            ElemMat::matrix_view mat_icomp = elmat.block( iblock + icoor, jblock + jcoor );
            // Loop on rows
            for ( i = 0; i < fe1.nbFEDof(); ++i )
            {
                // Loop on columns
                for ( j = 0; j < fe2.nbFEDof(); ++j )
                {
                    sum = 0.0;
                    for ( ig = 0; ig < bdfe.nbQuadPt(); ig++ )
                        sum += phid1[ i ][ icoor ][ ig ] * phid2[ j ][ jcoor ][ ig ] * bdfe.weightMeas( ig );
                    mat_icomp( i, j ) += coef * sum;
                }
            }
        }
    }

}

void ipstab_bagrad( const Real coef, ElemMat& elmat,
                    const CurrentFE& fe1, const CurrentFE& fe2,
                    const ElemVec& beta, const CurrentBdFE& bdfe,
                    int iblock, int jblock )
{

    // Interior penalty stabilization:
    // coef < |\beta . n|^2 / |\beta| \grad p1, \grad q2 >

    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    Real sum, sum1, sum2;
    UInt icoor,jcoor;
    UInt i, j, ig;
    Real phid1[ fe1.nbFEDof() ][ fe1.nbCoor() ][ bdfe.nbQuadPt() ];
    Real phid2[ fe2.nbFEDof() ][ fe2.nbCoor() ][ bdfe.nbQuadPt() ];

    std::vector<Real> x(3), rx1(3), drp1(3), rx2(3), drp2(3);
    std::vector<Real> b1(3), b2(3);

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    //
    // convection velocity term |\beta . n|^2 / |\beta|
    // on the boundary quadrature points
    //
    Real ba2[ bdfe.nbQuadPt() ];

    for ( ig = 0; ig < bdfe.nbQuadPt(); ig++ )
    {
        sum1 = 0;
        sum2 = 0;
        for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
        {
            for ( i = 0; i < bdfe.nbCoor(); ++i )
            {
                Real betaLoc = bdfe.phi( i, ig ) *
                    beta.vec() [ icoor * bdfe.nbCoor() + i ];
                sum1 += betaLoc * bdfe.normal(icoor, ig);
                sum2 += betaLoc * betaLoc;
            }
        }
        ba2[ ig ] = sum2 == 0 ? 0 : sum1 * sum1 / pow( sum2, 0.5 );
    }

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac( icoor, jcoor, 0 ) * drp2[ jcoor ];
                }
                phid1[ i ][ icoor ][ ig ] = sum1;
                phid2[ i ][ icoor ][ ig ] = sum2;
            }
        }
    }



    // Loop on rows
    for ( i = 0; i < fe1.nbFEDof(); ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbFEDof(); ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt() ; ++ig )
                    sum += ba2[ ig ] *
                           phid1[ i ][ icoor ][ ig ] *
                           phid2[ j ][ icoor ][ ig ] *
                           bdfe.weightMeas( ig );
            mat( i, j ) += coef * sum;
        }
    }

}


// coef < |\beta . n| \grad p1, \grad q2 >
// p1   lives in fe1
// q2   lives in fe2
// beta lives in fe3


void ipstab_bagrad( const Real         coef,
                    ElemMat&           elmat,
                    const CurrentFE&   fe1,
                    const CurrentFE&   fe2,
                    const CurrentFE&   fe3,
                    const ElemVec&     beta,
                    const CurrentBdFE& bdfe,
                    int iblock, int    jblock )
{

    // Interior penalty stabilization:
    // coef < |\beta.n| \grad p1, \grad q2 >

    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    Real sum, sum1, sum2;
    UInt icoor,jcoor;
    UInt i, j, ig;
    Real phid1[ fe1.nbFEDof() ][ fe1.nbCoor() ][ bdfe.nbQuadPt() ];
    Real phid2[ fe2.nbFEDof() ][ fe2.nbCoor() ][ bdfe.nbQuadPt() ];

    std::vector<Real> x(3), rx1(3), drp1(3), rx2(3), drp2(3);
    std::vector<Real> b1(3), b2(3);

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    //
    // convection velocity term |\beta . n|
    // on the boundary quadrature points
    //
    Real bn[ bdfe.nbQuadPt() ];

    for ( ig = 0; ig < bdfe.nbQuadPt(); ig++ )
    {
        sum1 = 0;
        sum2 = 0;

        for ( icoor = 0; icoor < fe3.nbCoor(); ++icoor )
        {
            for ( i = 0; i < fe3.nbFEDof(); ++i )
            {
                Real betaLoc = fe3.phi( i, ig ) *
                               beta.vec() [ icoor*fe3.nbFEDof() + i ];
                sum1 += betaLoc * bdfe.normal(icoor, ig);
            }
        }
        bn[ ig ] = std::abs(sum1);
    }

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac( icoor, jcoor, 0 ) * drp2[ jcoor ];
                }
                phid1[ i ][ icoor ][ ig ] = sum1;
                phid2[ i ][ icoor ][ ig ] = sum2;
            }
        }
    }



    // Loop on rows
    for ( i = 0; i < fe1.nbFEDof(); ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbFEDof(); ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor(); ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt() ; ++ig )
                    sum += bn[ ig ] *
                           phid1[ i ][ icoor ][ ig ] *
                           phid2[ j ][ icoor ][ ig ] *
                           bdfe.weightMeas( ig );
            mat( i, j ) += coef * sum;
        }
    }

}


void stiff( Real coef, ElemMat& elmat, const CurrentFE& fe,
            int iblock, int jblock )
/*
  Stiffness matrix: coef*\int grad v_i . grad v_j
*/
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt iloc, jloc;
    UInt i, icoor, ig;
    double s, coef_s;
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig )
                     * fe.weightDet( ig );
        }
        mat( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                     fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat( iloc, jloc ) += coef_s;
        mat( jloc, iloc ) += coef_s;
    }
}


void stiff( Real coef, Real ( *fct ) ( Real, Real, Real ), ElemMat& elmat,
            const CurrentFE& fe, int iblock, int jblock )
/*
  Stiffness matrix: coef*\int grad v_i . grad v_j
*/
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt iloc, jloc;
    UInt i, icoor, ig;
    double s, coef_s, coef_f;
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            coef_f = fct( fe.quadPt( ig, 0 ), fe.quadPt( ig, 1 ), fe.quadPt( ig, 2 ) );
            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                s += coef_f * fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig )
                     * fe.weightDet( ig );
        }
        mat( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            coef_f = fct( fe.quadPt( ig, 0 ), fe.quadPt( ig, 1 ), fe.quadPt( ig, 2 ) );
            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                s += coef_f * fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                     fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat( iloc, jloc ) += coef_s;
        mat( jloc, iloc ) += coef_s;
    }
}
//

void stiff( Real coef, ElemMat& elmat, const CurrentFE& fe,
            int iblock, int jblock, int nb )
/*
  Stiffness matrix: coef*\int grad v_i . grad v_j (nb blocks on the diagonal, nb>1)
*/
{


    Matrix mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );

    UInt iloc, jloc;
    UInt i, icoor, ig;
    double s, coef_s;
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig )
                     * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                     fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp( iloc, jloc ) += coef_s;
        mat_tmp( jloc, iloc ) += coef_s;
    }
    // copy on the components
    for ( int icomp = 0; icomp < nb; icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst( i );
            jloc = fe.patternSecond( i );
            mat_icomp( iloc, jloc ) += mat_tmp( iloc, jloc );
            mat_icomp( jloc, iloc ) += mat_tmp( jloc, iloc );
        }
    }
}


void stiff( const std::vector<Real>& coef, ElemMat& elmat, const CurrentFE& fe,
            int iblock, int jblock, int nb )
/*
  Stiffness matrix: coef*\int grad v_i . grad v_j (nb blocks on the diagonal, nb>1)
  with coef given in a vector (one element per quadrature point)
*/
{

    Matrix mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );

    UInt iloc, jloc;
    UInt i, icoor, ig;
    double s;//, coef_s;
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                s += coef[ig] *fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig )
                     * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) += s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                s += coef[ig] * fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                     fe.weightDet( ig );
        }
        mat_tmp( iloc, jloc ) += s;
        mat_tmp( jloc, iloc ) += s;
    }
    // copy on the components
    for ( int icomp = 0; icomp < nb; icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst( i );
            jloc = fe.patternSecond( i );
            mat_icomp( iloc, jloc ) += mat_tmp( iloc, jloc );
            mat_icomp( jloc, iloc ) += mat_tmp( jloc, iloc );
        }
    }
}


// VC - December 2004
//
void stiff_curl( Real coef, ElemMat& elmat, const CurrentFE& fe,
                 int iblock, int jblock, int /*nb*/ )


/*
  Stiffness matrix: coef*\int curl v_i . curl v_j (nb blocks on the diagonal, nb>1)
*/
{


    Matrix mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp11( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp11 = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp12( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp12 = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp13( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp13 = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp21( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp21 = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp22( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp22 = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp23( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp23 = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp31( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp31 = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp32( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp32 = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp33( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp33 = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );



    UInt iloc, jloc;
    UInt i, ig;
    double s, coef_s;

    // diagonal 11
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer( iloc, 1, ig ) * fe.phiDer( iloc, 1, ig ) * fe.weightDet( ig )
                + fe.phiDer( iloc, 2, ig ) * fe.phiDer( iloc, 2, ig ) * fe.weightDet( ig ) ;
        }
        mat_tmp11( iloc, iloc ) += coef * s;
    }
    // extra diagonal 11
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer( iloc, 1, ig ) * fe.phiDer( jloc, 1, ig ) * fe.weightDet( ig )
                + fe.phiDer( iloc, 2, ig ) * fe.phiDer( jloc, 2, ig ) * fe.weightDet( ig )       ;
        }
        coef_s = coef * s;
        mat_tmp11( iloc, jloc ) += coef_s;
        mat_tmp11( jloc, iloc ) += coef_s;
    }

    // diagonal 12
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 1, ig ) * fe.phiDer( iloc, 0, ig ) * fe.weightDet( ig );
        }
        mat_tmp12( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 12
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 1, ig ) * fe.phiDer( jloc, 0, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp12( iloc, jloc ) -= coef_s;
        mat_tmp12( jloc, iloc ) -= coef_s;
    }

    // diagonal 13
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 2, ig ) * fe.phiDer( iloc, 0, ig ) * fe.weightDet( ig );
        }
        mat_tmp13( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 13
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 2, ig ) * fe.phiDer( jloc, 0, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp13( iloc, jloc ) -= coef_s;
        mat_tmp13( jloc, iloc ) -= coef_s;
    }

    // diagonal 21
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 0, ig ) * fe.phiDer( iloc, 1, ig ) * fe.weightDet( ig );
        }
        mat_tmp21( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 21
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 0, ig ) * fe.phiDer( jloc, 1, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp21( iloc, jloc ) -= coef_s;
        mat_tmp21( jloc, iloc ) -= coef_s;
    }

    // diagonal 22
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer( iloc, 0, ig ) * fe.phiDer( iloc, 0, ig ) * fe.weightDet( ig )
                + fe.phiDer( iloc, 2, ig ) * fe.phiDer( iloc, 2, ig ) * fe.weightDet( ig ) ;
        }
        mat_tmp22( iloc, iloc ) += coef * s;
    }
    // extra diagonal 22
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer( iloc, 0, ig ) * fe.phiDer( jloc, 0, ig ) * fe.weightDet( ig )
                + fe.phiDer( iloc, 2, ig ) * fe.phiDer( jloc, 2, ig ) * fe.weightDet( ig )       ;
        }
        coef_s = coef * s;
        mat_tmp22( iloc, jloc ) += coef_s;
        mat_tmp22( jloc, iloc ) += coef_s;
    }

    // diagonal 23
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 2, ig ) * fe.phiDer( iloc, 1, ig ) * fe.weightDet( ig );
        }
        mat_tmp23( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 23
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 2, ig ) * fe.phiDer( jloc, 1, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp23( iloc, jloc ) -= coef_s;
        mat_tmp23( jloc, iloc ) -= coef_s;
    }

    // diagonal 31
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 0, ig ) * fe.phiDer( iloc, 2, ig ) * fe.weightDet( ig );
        }
        mat_tmp31( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 31
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 0, ig ) * fe.phiDer( jloc, 2, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp31( iloc, jloc ) -= coef_s;
        mat_tmp31( jloc, iloc ) -= coef_s;
    }

    // diagonal 32
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 1, ig ) * fe.phiDer( iloc, 2, ig ) * fe.weightDet( ig );
        }
        mat_tmp32( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 32
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer( iloc, 1, ig ) * fe.phiDer( jloc, 2, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp32( iloc, jloc ) -= coef_s;
        mat_tmp32( jloc, iloc ) -= coef_s;
    }

    // diagonal 33
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer( iloc, 0, ig ) * fe.phiDer( iloc, 0, ig ) * fe.weightDet( ig )
                + fe.phiDer( iloc, 1, ig ) * fe.phiDer( iloc, 1, ig ) * fe.weightDet( ig ) ;
        }
        mat_tmp33( iloc, iloc ) += coef * s;
    }
    // extra diagonal 33
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer( iloc, 1, ig ) * fe.phiDer( jloc, 1, ig ) * fe.weightDet( ig )
                + fe.phiDer( iloc, 2, ig ) * fe.phiDer( jloc, 2, ig ) * fe.weightDet( ig )       ;
        }
        coef_s = coef * s;
        mat_tmp33( iloc, jloc ) += coef_s;
        mat_tmp33( jloc, iloc ) += coef_s;
    }

    ElemMat::matrix_view mat_icomp = elmat.block( iblock + 0, jblock + 0 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) += mat_tmp11( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) += mat_tmp11( iloc, jloc );
        mat_icomp( jloc, iloc ) += mat_tmp11( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 0, jblock + 1 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp12( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp12( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp12( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 0, jblock + 2 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp13( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp13( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp13( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 1, jblock + 0 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp21( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp21( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp21( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 1, jblock + 1 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) += mat_tmp22( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) += mat_tmp22( iloc, jloc );
        mat_icomp( jloc, iloc ) += mat_tmp22( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 1, jblock + 2 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp23( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp23( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp23( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 2, jblock + 0 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp31( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp31( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp31( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 2, jblock + 1 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp32( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp32( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp32( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 2, jblock + 2 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) += mat_tmp33( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) += mat_tmp33( iloc, jloc );
        mat_icomp( jloc, iloc ) += mat_tmp( jloc, iloc );
    }
}


/*
  Stiffness matrix: coef * ( div u , div v )
*/
void stiff_div( Real coef, ElemMat& elmat, const CurrentFE& fe )
{

    double s;

    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

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



/*
  Stiffness matrix: coef * ( [\grad u^k]^T \grad d : \grad v  )
*/
void stiff_dergradbis( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

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




/*
  Stiffness matrix: coef * ( [\grad u]^T \grad u^k [\grad u^k]^T \grad u : \grad v  ) for Newton on St-Venant
*/
void stiff_dergrad( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

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





//
// coef * ( \tr { [\grad u^k]^T \grad u }, \div v  ) for Newton on St-Venant
//
//
void stiff_derdiv( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe )
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

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

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




void stiff_strain( Real coef, ElemMat& elmat, const CurrentFE& fe )
/*
  Stiffness matrix: coef * ( e(u) , e(v) )
*/
{
    double s;
    double tmp = coef * 0.5;

    ElemMat::matrix_type mat_tmp( fe.nbFEDof(), fe.nbFEDof() );

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
        ElemMat::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
    }

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {
            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );
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



void mass_divw( Real coef, const ElemVec& w_loc, ElemMat& elmat, const CurrentFE& fe,
                int iblock, int jblock, UInt nb )
/*
  modified mass matrix: ( div w u,v )
*/
{

    Matrix mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );

    UInt i, icomp, ig, icoor, iloc, jloc;
    Real s, coef_s, divw[ fe.nbQuadPt() ];

    // divw at quadrature nodes
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        divw[ ig ] = 0.0;
        for ( icoor = 0; icoor < fe.nbCoor(); ++icoor )
            for ( i = 0; i < fe.nbFEDof(); ++i )
                divw[ ig ] += fe.phiDer( i, icoor, ig ) * w_loc.vec() [ i + icoor * fe.nbFEDof() ];
    }

    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += divw[ ig ] * fe.phi( iloc, ig ) * fe.phi( iloc, ig ) * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
            s += divw[ ig ] * fe.phi( iloc, ig ) * fe.phi( jloc, ig ) * fe.weightDet( ig );
        coef_s = coef * s;
        mat_tmp( iloc, jloc ) += coef_s;
        mat_tmp( jloc, iloc ) += coef_s;
    }
    // copy on the components
    for ( icomp = 0; icomp < nb; icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst( i );
            jloc = fe.patternSecond( i );
            mat_icomp( iloc, jloc ) += mat_tmp( iloc, jloc );
            mat_icomp( jloc, iloc ) += mat_tmp( jloc, iloc );
        }
    }
}


void mass_divw(const std::vector<Real>& coef, const ElemVec& w_loc, ElemMat& elmat, const CurrentFE& fe,
               int iblock, int jblock, UInt nb )
/*
  modified mass matrix: ( div w u,v )
*/
{

    Matrix mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );

    UInt i, icomp, ig, icoor, iloc, jloc;
    Real s, divw[ fe.nbQuadPt() ]; // , coef_s

    // divw at quadrature nodes
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        divw[ ig ] = 0.0;
        for ( icoor = 0; icoor < fe.nbCoor(); ++icoor )
            for ( i = 0; i < fe.nbFEDof(); ++i )
                divw[ ig ] += fe.phiDer( i, icoor, ig ) * w_loc.vec() [ i + icoor * fe.nbFEDof() ];
    }

    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += coef[ig] * divw[ ig ] * fe.phi( iloc, ig ) * fe.phi( iloc, ig ) * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) += s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
            s += coef[ig]* divw[ ig ] * fe.phi( iloc, ig ) * fe.phi( jloc, ig ) * fe.weightDet( ig );
        mat_tmp( iloc, jloc ) += s;
        mat_tmp( jloc, iloc ) += s;
    }
    // copy on the components
    for ( icomp = 0; icomp < nb; icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst( i );
            jloc = fe.patternSecond( i );
            mat_icomp( iloc, jloc ) += mat_tmp( iloc, jloc );
            mat_icomp( jloc, iloc ) += mat_tmp( jloc, iloc );
        }
    }
}




void mass_gradu( Real coef, const ElemVec& u0_loc, ElemMat& elmat, const CurrentFE& fe )
/*
  modified mass matrix: ( grad u0 u,v )
*/
{

    UInt ig, icoor, jcoor, i, j;
    Real s;
    Real gu0[ fe.nbQuadPt() ][ fe.nbCoor() ][ fe.nbCoor() ];


    //
    // grad u0 at quadrature nodes
    //
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        for ( icoor = 0; icoor < fe.nbCoor(); ++icoor )
            for ( jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
            {
                gu0[ ig ][ icoor ][ jcoor ] = 0.0;
                for ( i = 0; i < fe.nbFEDof(); ++i )
                    gu0[ ig ][ icoor ][ jcoor ] += fe.phiDer( i, jcoor, ig ) * u0_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        for ( jcoor = 0; jcoor < fe.nbCoor(); ++jcoor )
        {

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

            for ( i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( ig = 0; ig < fe.nbQuadPt(); ++ig )
                        s += gu0[ ig ][ icoor ][ jcoor ] * fe.phi( i, ig ) * fe.phi( j, ig ) * fe.weightDet( ig );
                    mat( i, j ) += coef * s;
                }
            }
        }
    }
}


//
//
// \! Streamline diffusion
//
//
void stiff_sd( Real coef, const ElemVec& vec_loc, ElemMat& elmat, const CurrentFE& fe, const CurrentFE& fe2,
               int iblock, int jblock, int nb )
/*
  Stiffness matrix for SD Stabilization
*/
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt iloc, jloc;
    UInt i, icoor, ig, jcoor;
    double s, coef_s, coef_v[ nDimensions ];
    //    int nbN1=fe.nbFEDof();
    UInt nbN2 = fe2.nbFEDof();
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < nDimensions; icoor++ )
                coef_v[ icoor ] = 0.;

            // computation of the convection term in the quadrature nodes
            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
            {
                for ( UInt iloc = 0; iloc < nbN2; iloc++ )
                    coef_v[ icoor ] += vec_loc.vec() [ iloc + icoor * nbN2 ] * fe2.phi( iloc, ig );
            }

            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
            {
                for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                {
                    s += coef_v[ icoor ] * fe.phiDer( iloc, icoor, ig ) * coef_v[ jcoor ] * fe.phiDer( iloc, jcoor, ig )
                         * fe.weightDet( ig );
                }
            }
        }
        mat( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < nDimensions; icoor++ )
                coef_v[ icoor ] = 0.;

            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
            {
                for ( UInt iloc = 0; iloc < nbN2; iloc++ )
                    coef_v[ icoor ] += vec_loc.vec() [ iloc + icoor * nbN2 ] * fe2.phi( iloc, ig );
            }

            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
            {
                for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                {
                    s += coef_v[ icoor ] * fe.phiDer( iloc, icoor, ig ) * coef_v[ jcoor ] * fe.phiDer( jloc, jcoor, ig )
                         * fe.weightDet( ig );
                }
            }
        }
        coef_s = coef * s;
        mat( iloc, jloc ) += coef_s;
        mat( jloc, iloc ) += coef_s;
    }
    // copy on the other components (if necessary, i.e. if nb>1)
    for ( int icomp = 1; icomp < nb; icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst( i );
            jloc = fe.patternSecond( i );
            mat_icomp( iloc, jloc ) += mat( iloc, jloc );
            mat_icomp( jloc, iloc ) += mat( jloc, iloc );
        }
    }
}


//
void grad( const int icoor, Real coef, ElemMat& elmat,
           const CurrentFE& fe_u, const CurrentFE& fe_p,
           int iblock, int jblock )
/*
  \int q_j \frac{\partial v_i}{\partial x_icoor}
*/
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt ig;
    UInt i, j;
    double s;
    for ( i = 0; i < fe_u.nbFEDof(); i++ )
    {
        for ( j = 0; j < fe_p.nbFEDof(); j++ )
        {
            s = 0;
            for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
                // Be careful the minus is here and not in coef!!!!

                // wrong for different quadrules of fe_u and fe_p !!   Martin P.
                // s -= fe_p.phi(j,ig)*fe_u.phiDer(i,icoor,ig)*fe_u.weightDet(ig);

                s -= fe_p.refFE().phi( j, fe_u.quadRule().quadPointCoor( ig, 0 ), fe_u.quadRule().quadPointCoor( ig, 1 ),
                                       fe_u.quadRule().quadPointCoor( ig, 2 ) ) * fe_u.phiDer( i, icoor, ig ) * fe_u.weightDet( ig );
            mat( i, j ) += coef * s;
        }
    }
}

void div( const int icoor, Real coef, ElemMat& elmat,
          const CurrentFE& fe_u, const CurrentFE& fe_p,
          int iblock, int jblock )
/*
  \int q_i \frac{\partial v_j}{\partial x_icoor}
*/
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt ig;
    UInt i, j;
    double s;
    for (i = 0; i < fe_u.nbFEDof(); i++)
    {
        for (j = 0; j < fe_p.nbFEDof(); j++)
        {
            s = 0;
            for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
                s -= fe_u.phi(i, ig) * fe_p.phiDer(j, icoor, ig) * fe_u.weightDet( ig );
            mat( i, j ) += coef * s;
        }
    }
}

void grad_div( Real coef_grad, Real coef_div, ElemMat& elmat,
               const CurrentFE& fe_u, const CurrentFE& fe_p,
               int block_pres )
/*
  \int q_j \frac{\partial v_i}{\partial x_icoor}
*/
{
    double s;
    int iblock = block_pres - nDimensions;
    for ( UInt icoor = 0; icoor < 3; icoor++ )
    {
        ElemMat::matrix_view mat_grad = elmat.block( iblock + icoor, block_pres );
        ElemMat::matrix_view mat_div = elmat.block( block_pres , iblock + icoor );
        for ( UInt i = 0; i < fe_u.nbFEDof(); i++ )
        {
            for ( UInt j = 0; j < fe_p.nbFEDof(); j++ )
            {
                s = 0;
                for ( UInt ig = 0; ig < fe_u.nbQuadPt(); ig++ )
                    s -= fe_p.phi( j, ig ) * fe_u.phiDer( i, icoor, ig ) * fe_u.weightDet( ig );
                mat_grad( i, j ) += coef_grad * s;
                mat_div( j, i ) += coef_div * s;
            }
        }
    }
}
//
void stab_stokes( Real visc, Real coef_stab, ElemMat& elmat,
                  const CurrentFE& fe, int block_pres )
{
    ElemMat::matrix_view mat = elmat.block( block_pres, block_pres );
    Real s, h = fe.diameter();
    Real fh2 = coef_stab * h * h / ( 2 * visc );
    for ( UInt i = 0; i < fe.nbFEDof(); i++ )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); j++ )
        {
            s = 0;
            for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
            {
                for ( UInt icoor = 0; icoor < 3; icoor++ )
                {
                    s += fe.phiDer( i, icoor, ig ) * fe.phiDer( j, icoor, ig ) * fe.weightDet( ig );
                }
            }
            mat( i, j ) -= fh2 * s;
        }
    }
}

/*
 * Fixed by Umberto Villa,  Jan 2010
 */
void advection( Real coef, ElemVec& vel,
                ElemMat& elmat, const CurrentFE& fe, int iblock, int jblock, int nb )
{
    Matrix mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    Real v_grad, s;
    Matrix v( fe.nbQuadPt(), nDimensions );

    //Evaluate the advective field at the quadrature nodes
    for ( UInt icoor = 0; icoor < nDimensions; icoor++ )
    {
        ElemVec::vector_view velicoor = vel.block( icoor );
        for ( UInt iq = 0; iq < fe.nbQuadPt(); iq++ )
        {
            s = 0.;
            for ( UInt k = 0; k < fe.nbFEDof(); k++ )
                s += velicoor( k ) * fe.phi( k, iq ); // velocity on the intgt point
            v(iq, icoor) = s;
        }
    }

    //Assemble the local matrix
    for ( UInt i = 0; i < fe.nbFEDof(); i++ )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); j++ )
        {
            s = 0.;
            for ( UInt iq = 0; iq < fe.nbQuadPt(); iq++ )
            {
                v_grad = 0.;
                for ( int icoor = 0; icoor < ( int ) nDimensions; icoor++ )
                    v_grad += v(iq, icoor) * fe.phiDer( j, icoor, iq );

                s += v_grad * fe.phi( i, iq ) * fe.weightDet( iq );
            }
            mat_tmp( i, j ) = s*coef;
        }
    }

    // copy on the components
    for ( int icomp = 0; icomp < nb; icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( UInt i = 0; i < fe.nbDiag(); i++ )
        {
            for ( UInt j = 0; j < fe.nbDiag(); j++ )
            {
                mat_icomp( i, j ) += mat_tmp( i, j );
            }
        }
    }
}


void grad( const int icoor, const ElemVec& vec_loc, ElemMat& elmat,
           const CurrentFE& fe1, const CurrentFE& fe2,
           int iblock, int jblock )
/*
  \int q_j \frac{\partial v_i}{\partial x_icoor}
*/
{
    //
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    if ( iblock == jblock )
    {
        UInt iq;
        UInt i, j;
        double s, coef;
        UInt nbN1 = fe1.nbFEDof();
        for ( i = 0; i < nbN1; i++ )
        {
            for ( j = 0; j < fe2.nbFEDof(); j++ )
            {
                s = 0;
                for ( iq = 0; iq < fe1.nbQuadPt(); iq++ )
                {
                    coef = 0;
                    for ( UInt iloc = 0; iloc < nbN1; iloc++ )
                        coef += vec_loc.vec() [ iloc + icoor * nbN1 ] * fe1.phi( iloc, iq );

                    s += coef * fe2.phi( i, iq ) * fe1.phiDer( j, icoor, iq ) * fe1.weightDet( iq );
                } // Loop on quadrature nodes

                mat( i, j ) += s;
            } //Loop on j
        } // Loop on i
    } // if

}

//
// \! Gradient operator in the skew-symmetric form for NS Problems
//  A. Veneziani - December 2002
// \!
void grad_ss( const int icoor, const ElemVec& vec_loc, ElemMat& elmat,
              const CurrentFE& fe1, const CurrentFE& fe2,
              int iblock, int jblock )
/*
  \int vloc(icoor) \frac{\partial v_i}{\partial x_icoor} v_j + 1/2*\frac{\partial v_icoor}{\partial x_icoor} v_i v_j
*/
{
    //
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    if ( iblock == jblock )
    {
        UInt iq;
        UInt i, j;
        double s, coef, coef_div;
        UInt nbN1 = fe1.nbFEDof();
        for ( i = 0; i < nbN1; i++ )
        {
            for ( j = 0; j < fe2.nbFEDof(); j++ )
            {
                s = 0;
                for ( iq = 0; iq < fe1.nbQuadPt(); iq++ )
                {
                    coef = 0;
                    coef_div = 0;

                    for ( UInt iloc = 0; iloc < nbN1; iloc++ )
                    {
                        coef += vec_loc.vec() [ iloc + icoor * nbN1 ] * fe1.phi( iloc, iq );
                        coef_div += vec_loc.vec() [ iloc + icoor * nbN1 ] * fe1.phiDer( iloc, icoor, iq );
                    }

                    s += ( coef * fe1.phiDer( j, icoor, iq ) + 0.5 * coef_div * fe1.phi( j, iq ) ) * fe2.phi( i, iq ) * fe1.weightDet( iq );
                } // Loop on quadrature nodes

                mat( i, j ) += s;
            } //Loop on j
        } // Loop on i
    } // if

}

// /!
// Gradient operator where the convective term is based on a local vector
// living on the basis given by fe3
// It is useful for advection diffusion problems driven by a NS problem
// !/
void grad( const int icoor,
           const ElemVec& vec_loc,
           ElemMat& elmat,
           const CurrentFE& fe1,
           const CurrentFE& fe2,
           const CurrentFE& fe3,
           int iblock, int jblock )
/*
  \int q_j \frac{\partial v_i}{\partial x_icoor}
*/
{
    //
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    if ( iblock == jblock )
    {
        UInt iq;
        UInt i, j;
        double s, coef;
        UInt nbN1 = fe1.nbFEDof();
        UInt nbN3 = fe3.nbFEDof();

        for ( i = 0; i < nbN1; i++ )
        {
            for ( j = 0; j < fe2.nbFEDof(); j++ )
            {
                s = 0;
                for ( iq = 0; iq < fe1.nbQuadPt(); iq++ )
                {
                    coef = 0;

                    for ( UInt iloc = 0; iloc < nbN3; iloc++ )
                        coef += vec_loc.vec() [iloc + icoor*nbN3] * fe3.phi(iloc, iq);

                    s += coef * fe2.phi(i, iq) * fe1.phiDer(j, icoor, iq) * fe1.weightDet(iq);
                } // Loop on quadrature nodes

                mat(i, j) += s;
            } //Loop on j
        } // Loop on i
    } // if

}


// Convective term with the velocity given in the quadrature nodes
void grad( const int& icoor,
           const std::vector<Real>& localVector,
           ElemMat& elmat,
           const CurrentFE& currentFE1,
           const CurrentFE& currentFE2,
           const int& iblock,
           const int& jblock)
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    // This term concerns only the diagonal blocks (same components)
    if ( iblock == jblock )
    {

        for ( UInt iNode1(0); iNode1 < currentFE1.nbFEDof(); iNode1++ )
        {
            for ( UInt jNode2(0); jNode2 < currentFE2.nbFEDof(); jNode2++ )
            {
                Real sum(0.0);
                for ( UInt iq(0); iq < currentFE1.nbQuadPt(); iq++ )
                {
                    // Velocity in the quadrature node, component icoor
                    double coef(localVector[icoor*currentFE1.nbQuadPt() + iq]);

                    sum += coef
                           * currentFE2.phi(iNode1, iq)
                           * currentFE1.phiDer(jNode2, icoor, iq)
                           * currentFE1.weightDet(iq);
                } // Loop on quadrature nodes

                mat(iNode1, jNode2) += sum;
            } //Loop on j
        } // Loop on i
    } // if

}
//
//
/*
//////////////////////
// NONLINEAR TERMS
/////////////////////
//
//
//--------------------
// Jacobian: 2*\Sum V_k \Int \phi_k \phi_i \phi_j
//    and
// Vector F(V) = \Sum V_k \Sum V_j \Int \phi_k \phi_i \phi_j
//-------------------
void quad(std::vector<Real> coef, ElemMat& elmat, ElemVec& elvec,
const CurrentFE& fe,int iblock=0,int jblock=0)
{
ElemMat::matrix_view mat = elmat.block(iblock,jblock);
int i,ig,iq,siz;
int iloc,jloc,qloc;
Real s,coef_s;
siz=coef.size();
ASSERT(siz==fe.nbDiag,
"Error in building Local Matrix of the quadratic term");
//
// diagonal
//
for(i=0;i<fe.nbDiag();i++){
iloc = fe.patternFirst(i);
s = 0;
for (iq=0;i<siz;++i){
qloc = fe.patternFirst(iq);
for(ig=0;ig<fe.nbQuadPt();ig++)
s += coef(iq)*fe.phi(qloc,ig)*fe.phi(iloc,ig)*fe.phi(iloc,ig)*
fe.weightDet(ig);
}
mat(iloc,iloc) += 2*s;
elvec(iloc) += coef(iloc)*s;
}
//
// extra diagonal
//
for(i=fe.nbDiag();i<fe.nbDiag()+fe.nbUpper();i++){
iloc = fe.patternFirst(i);
jloc = fe.patternSecond(i);
s = 0;
for (iq=0;i<siz;++i){
qloc = fe.patternFirst(iq);
for(ig=0;ig<fe.nbQuadPt();ig++)
s += coef(iq)*fe.phi(qloc,ig)*
fe.phi(iloc,ig)*fe.phi(jloc,ig)*fe.weightDet(ig);
}
mat(iloc,jloc) += 2*s;
mat(jloc,iloc) += 2*s;
}
}
*/
//----------------------------------------------------------------------
//                      Element vector operator
//----------------------------------------------------------------------
void source( Real constant, ElemVec& elvec, const CurrentFE& fe, int iblock )
{
    UInt i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;
    for ( i = 0; i < fe.nbFEDof(); i++ )
    {
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
            s += fe.phi( i, ig ) * fe.weightDet( ig );
        vec( i ) += constant * s;
    }
}

void source( Real coef, ElemVec& f, ElemVec& elvec, const CurrentFE& fe,
             int fblock, int eblock )
/*
  compute \int f \phi_i
  where f is given on the dof of this element
  fblock is the block of the values read in f
  eblock is the block where the result is writen in the elvec
*/
{
    UInt i, ig;
    ElemVec::vector_view vec = elvec.block( eblock );
    ElemVec::vector_view vecf = f.block( fblock );
    Real f_ig;

    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        f_ig = 0.;
        for ( i = 0; i < fe.nbFEDof(); i++ )
            f_ig += vecf( i ) * fe.phi( i, ig );
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {
            vec( i ) += coef * f_ig * fe.phi( i, ig ) * fe.weightDet( ig );
        }
    }
}


void source_mass(const std::vector<Real>& constant, ElemVec& elvec, const CurrentFE& currentFe, const int& iblock)
{
    ElemVec::vector_view vec = elvec.block( iblock );
    for (UInt iterNode(0); iterNode < currentFe.nbFEDof(); ++iterNode )
    {
        for ( UInt iterQuadNode(0); iterQuadNode < currentFe.nbQuadPt(); iterQuadNode++ )
        {
            vec(iterNode) += constant[iterQuadNode]
                             * currentFe.phi( iterNode, iterQuadNode )
                             * currentFe.weightDet( iterQuadNode );
        };
    };
};

void source_stiff(const std::vector<Real>& constant, ElemVec& elvec, const CurrentFE& currentFe, const int& iblock)
{
    ElemVec::vector_view vec = elvec.block( iblock );
    const UInt nbQuadPt(currentFe.nbQuadPt());
    for (UInt iterNode(0); iterNode < currentFe.nbFEDof(); ++iterNode )
    {
        for ( UInt iterQuadNode(0); iterQuadNode < nbQuadPt; iterQuadNode++ )
        {
            for (UInt iterGradComp(0); iterGradComp<currentFe.nbCoor(); ++iterGradComp)
            {
                vec(iterNode) += constant[iterQuadNode + iterGradComp*nbQuadPt]
                                 * currentFe.phiDer( iterNode,iterGradComp, iterQuadNode )
                                 * currentFe.weightDet( iterQuadNode );
            };
        };
    };
};



void source_divuq(Real alpha, ElemVec& uLoc,  ElemVec& elvec, const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock  )
{
    UInt i, j, ic, iq;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;

    for (i = 0; i < fe_p.nbFEDof(); i++)
    {
        s = 0;
        for (iq = 0; iq < fe_p.nbQuadPt(); ++iq )
            for (j = 0; j < fe_u.nbFEDof(); ++j)
                for (ic = 0; ic < nDimensions; ++ic)
                    s += uLoc[ic*fe_u.nbFEDof()+j]*fe_u.phiDer(j,ic,iq)*fe_p.phi(i,iq)* fe_p.weightDet( iq );

        vec( i ) += s*alpha;
    }
}


void source_gradpv(Real alpha, ElemVec& pLoc,  ElemVec& elvec, const CurrentFE& fe_p, const CurrentFE& fe_u, int iblock )
{
    UInt i, j, iq;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;

    for ( i = 0; i < fe_u.nbFEDof(); i++ )
    {
        s = 0;
        for (iq = 0; iq < fe_u.nbQuadPt(); ++iq )
            for (j = 0; j < fe_p.nbFEDof(); ++j)
                s += pLoc[j]*fe_p.phiDer(j,iblock,iq)*fe_u.phi(i,iq)*fe_u.weightDet( iq );
        vec( i ) += s*alpha;
    }
}




void source_fhn( Real coef_f, Real coef_a, ElemVec& u, ElemVec& elvec, const CurrentFE& fe,
                 int fblock, int eblock )
/*
  compute \int coef_f u(1-u)(u-coef_a) \phi_i
  (right-hand side for the Fitzhugh-Nagumo equations)
  where u is given on the dof of this element
  fblock is the block of the values read in f
  eblock is the block where the result is writen in the elvec
*/
{
    UInt i, ig;
    ElemVec::vector_view vec = elvec.block( eblock );
    ElemVec::vector_view vecu = u.block( fblock );
    Real f_ig;

    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        f_ig = 0.;
        for ( i = 0; i < fe.nbFEDof(); i++ )
            f_ig += vecu( i ) * fe.phi( i, ig );
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {
            vec( i ) += coef_f * f_ig * ( 1 - f_ig ) * ( f_ig - coef_a ) * fe.phi( i, ig ) * fe.weightDet( ig );
        }
    }

}

// coef * ( - \grad w^k :[I\div d - (\grad d)^T] u^k + convect^T[I\div d - (\grad d)^T] (\grad u^k)^T , v  ) for Newton FSI
//
// Remark: convect = u^n-w^k relative vel.
//
void source_mass1( Real coef, const ElemVec& uk_loc, const ElemVec& wk_loc, const ElemVec& convect_loc,
                   const ElemVec& d_loc, ElemVec& elvec, const CurrentFE& fe )
{
    Real B[ fe.nbCoor() ][ fe.nbCoor() ];                 // \grad (convect) at a quadrature point
    Real A[ fe.nbCoor() ][ fe.nbCoor() ];                 // I\div d - (\grad d)^T at a quadrature point
    Real aux[ fe.nbQuadPt() ];                        // grad (- w^k):[I\div d - (\grad d)^T] at  quadrature points
    Real uk[ fe.nbQuadPt() ][ fe.nbCoor() ];              // u^k quadrature points
    Real guk[ fe.nbQuadPt() ][ fe.nbCoor() ][ fe.nbCoor() ];  // \grad u^k at quadrature points
    Real convect[ fe.nbCoor() ];                      // convect at quadrature points
    Real convect_A[ fe.nbQuadPt() ][ fe.nbCoor() ];       // (convect)^T [I\div d - (\grad d)^T] at quadrature points

    Real s, sA, sB, sG;


    UInt icoor, jcoor, ig, i;
    // loop on quadrature points
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordindates
        for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // each compontent of uk at each quadrature points
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
                s += fe.phi( i, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
            uk[ ig ][ icoor ] = s;//uk_x(pt_ig), uk_y(pt_ig), uk_z(pt_ig)

            // each compontent of convect at this quadrature point
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
                s += fe.phi( i, ig ) * convect_loc.vec() [ i + icoor * fe.nbFEDof() ];
            convect[ icoor ] = s;


            // loop  on space coordindates
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                sB = 0.0;
                sA = 0.0;
                sG = 0.0;
                for ( i = 0; i < fe.nbFEDof(); i++ )
                {
                    sG += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at each quadrature point
                    sB -= fe.phiDer( i, jcoor, ig ) * wk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad (- w^k) at this quadrature point
                    sA -= fe.phiDer( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe.nbFEDof() ]; //  - (\grad d) ^T at this quadrature point
                }
                guk[ ig ][ icoor ][ jcoor ] = sG; // \grad u^k at each quadrature point
                B[ icoor ][ jcoor ] = sB; // \grad (convect) at this quadrature point
                A[ icoor ][ jcoor ] = sA; // -(\grad d) ^T at this quadrature point
            }
        }

        s = 0.0;
        for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            s -= A[ jcoor ][ jcoor ];  // \div d at this quadrature point ( - trace( A ) )

        for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            A[ jcoor ][ jcoor ] += s;  // I\div d - (\grad d)^T at this quadrature point (A+I(-tr(A)))

        s = 0;
        for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                s += B[ icoor ][ jcoor ] * A[ icoor ][ jcoor ]; // \grad (-w^k):[I\div d - (\grad d)^T] at each quadrature point
        aux[ ig ] = s;

        s = 0;
        for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
        {
            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                s += convect[ icoor ] * A[ icoor ][ jcoor ]; // convect^T [I\div d - (\grad d)^T]
            convect_A[ ig ][ jcoor ] = s;
        }
//         std::cout<<"aux "<<aux[ig]<<std::endl;

//     for(icoor=0; icoor<fe.nbCoor(); ++icoor)
//         for(jcoor=0; jcoor<fe.nbCoor(); ++jcoor)
//             {
//                 std::cout<<" Atrue ["<<icoor<<"]["<<jcoor<<"] = "<<A[icoor][jcoor]<<std::endl;
//             }

    }

    // At this point we have:
    //    v  \grad u^k at each quadrature point: guk
    //    v  convect^T [I\div d - (\grad d)^T] at each quadrature point: convect_A
    //    v  \grad (-w^k):[I\div d - (\grad d)^T]: aux


    //
    // Numerical integration
    //

    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
    {

        // the block iccor of the elementary vector
        ElemVec::vector_view vec = elvec.block( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
            {

                // \grad ( - w^k ):[I\div d - (\grad d)^T] \phi_i
                s += aux[ ig ] * uk[ ig ][ icoor ] * fe.phi( i, ig ) * fe.weightDet( ig );

                // convect^T [I\div d - (\grad d)^T] (\grad u^k)^T \phi_i
                for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                    s += convect_A[ ig ][ jcoor ] * guk[ ig ][ icoor ][ jcoor ] * fe.phi( i, ig ) * fe.weightDet( ig );
            }
            vec( i ) += coef * s;
        }
    }
}



//
// coef * ( \grad u^k dw, v  ) for Newton FSI
//
//
void source_mass2( Real coef, const ElemVec& uk_loc, const ElemVec& dw_loc,
                   ElemVec& elvec, const CurrentFE& fe )
{

    Real guk[ fe.nbCoor() ][ fe.nbCoor() ];      // \grad u^k at a quadrature point
    Real dw[ fe.nbCoor() ];                  // dw at a quadrature point
    Real aux[ fe.nbQuadPt() ][ fe.nbCoor() ];    // (\grad u^k)dw at each quadrature point
    Real s;

    UInt ig, icoor, jcoor, i;

    // loop on quadrature points
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // each compontent (icoor) of dw at this quadrature point
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
                s += fe.phi( i, ig ) * dw_loc.vec() [ i + icoor * fe.nbFEDof() ];
            dw[ icoor ] = s;

            // loop  on space coordinates
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                s = 0.0;
                for ( i = 0; i < fe.nbFEDof(); i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ] = s;
            }
        }

        // (\grad u^k)dw at each quadrature point
        for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {
            s = 0.0;
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                s += guk[ icoor ][ jcoor ] * dw[ jcoor ];
            aux[ ig ][ icoor ] = s;
        }
    }

    //
    // Numerical integration
    //

    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
    {

        ElemVec::vector_view vec = elvec.block( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
                s += aux[ ig ][ icoor ] * fe.phi( i, ig ) * fe.weightDet( ig );
            vec( i ) += coef * s;
        }
    }
}



//
// coef * ( \grad u^n :[2 I \div d - (\grad d)^T]  u^k , v  ) for Newton FSI
//
//
void source_mass3( Real coef, const ElemVec& un_loc, const ElemVec& uk_loc, const ElemVec& d_loc,
                   ElemVec& elvec, const CurrentFE& fe )
{
    Real B[ fe.nbCoor() ][ fe.nbCoor() ];                 // \grad u^n at a quadrature point
    Real A[ fe.nbCoor() ][ fe.nbCoor() ];                 // I\div d - (\grad d)^T at a quadrature point
    Real aux[ fe.nbQuadPt() ];                        //  \div d  \div u^n  + grad u^n:[I\div d - (\grad d)^T] at  quadrature points
    Real uk[ fe.nbQuadPt() ][ fe.nbCoor() ];              // u^k quadrature points

    Real s, sA, sB;


    UInt icoor, jcoor, ig, i;
    // loop on quadrature points
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordindates
        for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // each compontent of uk at each quadrature points
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
                s += fe.phi( i, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
            uk[ ig ][ icoor ] = s;


            // loop  on space coordindates
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                sB = 0.0;
                sA = 0.0;
                for ( i = 0; i < fe.nbFEDof(); i++ )
                {
                    sB += fe.phiDer( i, jcoor, ig ) * un_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^n at this quadrature point
                    sA -= fe.phiDer( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe.nbFEDof() ]; //  - (\grad d) ^T at this quadrature point
                }
                B[ icoor ][ jcoor ] = sB; // \grad u^n at this quadrature point
                A[ icoor ][ jcoor ] = sA; // -(\grad d) ^T at this quadrature point
            }
        }

        s = 0.0;
        for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            s -= A[ jcoor ][ jcoor ];  // \div d at this quadrature point ( - trace( A ) )

        for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            A[ jcoor ][ jcoor ] += 2 * s;  // 2 * I\div d - (\grad d)^T at this quadrature point

        s = 0;
        for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                s += B[ icoor ][ jcoor ] * A[ icoor ][ jcoor ]; // \grad u^n:[2 * I\div d - (\grad d)^T] at each quadrature point
        aux[ ig ] = s;
    }

    // At this point we have:
    //    v u^k at each quadrature point: uk
    //    v  \grad u^n:[ 2 * I\div d - (\grad d)^T]: aux


    //
    // Numerical integration
    //

    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
    {

        // the block iccor of the elementary vector
        ElemVec::vector_view vec = elvec.block( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {
            // loop on quadrature points
            s = 0;
            for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
                // \grad u^n:[2 * I\div d - (\grad d)^T] u^k \phi_i
                s += aux[ ig ] * uk[ ig ][ icoor ] * fe.phi( i, ig ) * fe.weightDet( ig );
            vec( i ) += coef * s;
        }
    }
}






//
// coef * ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  ) for Newton FSI
//
void source_stress( Real coef, Real mu, const ElemVec& uk_loc, const ElemVec& pk_loc,
                    const ElemVec& d_loc, ElemVec& elvec, const CurrentFE& fe_u,
                    const CurrentFE& fe_p )
{

    Real A[ fe_u.nbCoor() ][ fe_u.nbCoor() ];                 // I\div d - (\grad d)^T at a quadrature point
    Real guk[ fe_u.nbCoor() ][ fe_u.nbCoor() ];               // \grad u^k at a quadrature point
    Real sigma[ fe_u.nbCoor() ][ fe_u.nbCoor() ];             // [-p^k I + 2*mu e(u^k)] a quadrature point
    Real B[ fe_u.nbQuadPt() ][ fe_u.nbCoor() ][ fe_u.nbCoor() ];  // [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] at each quadrature point
    Real s, sA, sG, pk;

    UInt icoor, jcoor, kcoor, ig, i;

    // loop on quadrature points
    for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0; icoor < fe_u.nbCoor(); icoor++ )
        {

            // loop  on space coordindates
            for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
            {
                sA = 0.0;
                sG = 0.0;
                for ( i = 0; i < fe_u.nbFEDof(); i++ )
                {
                    sG += fe_u.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbFEDof() ]; //  \grad u^k at this quadrature point
                    sA -= fe_u.phiDer( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe_u.nbFEDof() ]; //  - (\grad d) ^T at this quadrature point
                }
                guk[ icoor ][ jcoor ] = sG;
                A[ icoor ][ jcoor ] = sA;
            }
        }

        s = 0.0;
        for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
            s -= A[ jcoor ][ jcoor ];  // \div d at a quadrature point ( - trace( A ) )

//         std::cout<<"div = "<< s <<std::endl;

        for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
            A[ jcoor ][ jcoor ] += s;  // I\div d  - (\grad d)^T

        pk = 0.0;
        for ( i = 0; i < fe_p.nbFEDof(); i++ )
            pk += fe_p.phi( i, ig ) * pk_loc.vec() [ i ]; // p^k at this quadrature point


        // sigma = [-p^k I + 2*mu e(u^k)] a quadrature point
        for ( icoor = 0; icoor < fe_u.nbCoor(); icoor++ )
        {
            for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
            {
                sigma[ icoor ][ jcoor ] = mu * ( guk[ icoor ][ jcoor ] + guk[ jcoor ][ icoor ] );
            }
            sigma[ icoor ][ icoor ] -= pk;
        }

        // [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] at each quadrature point
        for ( icoor = 0; icoor < fe_u.nbCoor(); icoor++ )
            for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
            {
                s = 0;
                for ( kcoor = 0; kcoor < fe_u.nbCoor(); kcoor++ )
                    s += sigma[ icoor ][ kcoor ] * A[ kcoor ][ jcoor ];
                B[ ig ][ icoor ][ jcoor ] = s;
            }
    }

//     for(icoor=0; icoor<fe_u.nbCoor(); ++icoor)
//         for(jcoor=0; jcoor<fe_u.nbCoor(); ++jcoor)
//             {
//                 std::cout<<" Atrue ["<<icoor<<"]["<<jcoor<<"] = "<<A[icoor][jcoor]<<std::endl;
//             }

//     for(icoor=0; icoor<fe_u.nbCoor(); ++icoor)
//         for(jcoor=0; jcoor<fe_u.nbCoor(); ++jcoor)
//             {
//                 double l=0.;
//                 for(int e=0; e<fe_u.nbQuadPt(); ++e)
//                     l+=B[e][icoor][jcoor];
//                 std::cout<<" Btrue ["<<icoor<<"]["<<jcoor<<"] = "<<l<<std::endl;
//             }

    //
    // Numerical integration
    //

    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0; icoor < fe_u.nbCoor(); icoor++ )
    {

        ElemVec::vector_view vec = elvec.block( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0; i < fe_u.nbFEDof(); i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )

                for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
                    s += B[ ig ][ icoor ][ jcoor ] * fe_u.phiDer( i, jcoor, ig ) * fe_u.weightDet( ig );
            vec( i ) += coef * s;
        }
    }
}



//
// + \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )
//
void source_stress2( Real coef, const ElemVec& uk_loc, const ElemVec& d_loc, ElemVec& elvec, const CurrentFE& fe_u )
{

    Real guk[ fe_u.nbCoor() ][ fe_u.nbCoor() ];               // \grad u^k at a quadrature point
    Real gd[ fe_u.nbCoor() ][ fe_u.nbCoor() ];                // \grad d at a quadrature point
    Real A[ fe_u.nbQuadPt() ][ fe_u.nbCoor() ][ fe_u.nbCoor() ];  // \grad u^k \grad d + [\grad d]^T[\grad u^k]^T  at each quadrature point
    Real su, sd, s;

    UInt icoor, jcoor, kcoor, ig, i;

    // loop on quadrature points
    for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0; icoor < fe_u.nbCoor(); icoor++ )
        {

            // loop  on space coordindates
            for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
            {
                su = 0.0;
                sd = 0.0;
                for ( i = 0; i < fe_u.nbFEDof(); i++ )
                {
                    su += fe_u.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbFEDof() ]; //  \grad u^k at this quadrature point
                    sd += fe_u.phiDer( i, jcoor, ig ) * d_loc.vec() [ i + icoor * fe_u.nbFEDof() ];  //  \grad d at this quadrature point
                }
                guk[ icoor ][ jcoor ] = su;
                gd[ icoor ][ jcoor ] = sd;
            }
        }


        for ( icoor = 0; icoor < fe_u.nbCoor(); icoor++ )
        {
            for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
            {
                s = 0;
                for ( kcoor = 0; kcoor < fe_u.nbCoor(); kcoor++ )
                    s += guk[ icoor ][ kcoor ] * gd[ kcoor ][ jcoor ] + gd[ kcoor ][ icoor ] * guk[ jcoor ][ kcoor ];
                A[ ig ][ icoor ][ jcoor ] = s;
            }
        }

//     for(icoor=0; icoor<fe_u.nbCoor(); ++icoor)
//         for(jcoor=0; jcoor<fe_u.nbCoor(); ++jcoor)
//             {
//                 std::cout<<" gdtrue ["<<icoor<<"]["<<jcoor<<"] = "<<gd[icoor][jcoor]<<std::endl;
//             }
//     for(icoor=0; icoor<fe_u.nbCoor(); ++icoor)
//         for(jcoor=0; jcoor<fe_u.nbCoor(); ++jcoor)
//             {
//                 std::cout<<" Atrue ["<<icoor<<"]["<<jcoor<<"] = "<<A[ig][icoor][jcoor]<<std::endl;
//             }
    }

    //
    // Numerical integration
    //
    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0; icoor < fe_u.nbCoor(); icoor++ )
    {

        ElemVec::vector_view vec = elvec.block( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0; i < fe_u.nbFEDof(); i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )

                for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
                    s += fe_u.phiDer( i, jcoor, ig ) * A[ ig ][ icoor ][ jcoor ] * fe_u.weightDet( ig );
            vec( i ) += coef * s;
        }
    }
}





//
// coef * (  (\grad u^k):[I\div d - (\grad d)^T] , q  ) for Newton FSI
//
void source_press( Real coef, const ElemVec& uk_loc, const ElemVec& d_loc, ElemVec& elvec,
                   const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock )
{

    Real A[ fe_u.nbCoor() ][ fe_u.nbCoor() ];     //  I\div d - (\grad d)^T at a quadrature point
    Real guk[ fe_u.nbCoor() ][ fe_u.nbCoor() ];   // \grad u^k at a quadrature point
    Real aux[ fe_u.nbQuadPt() ];              // grad u^k:[I\div d - (\grad d)^T] at each quadrature point
    ElemVec::vector_view vec = elvec.block( iblock );

    Real s, sA, sG;
    UInt icoor, jcoor, ig, i;


    // loop on quadrature points
    for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0; icoor < fe_u.nbCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
            {
                sA = 0.0;
                sG = 0.0;
                for ( i = 0; i < fe_u.nbFEDof(); i++ )
                {
                    sG += fe_u.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbFEDof() ]; //  \grad u^k at a quadrature point
                    sA -= fe_u.phiDer( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe_u.nbFEDof() ]; //  - (\grad d) ^T at a quadrature point
                }
                guk[ icoor ][ jcoor ] = sG;
                A[ icoor ][ jcoor ] = sA;
            }
        }

        s = 0.0;
        for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
            s -= A[ jcoor ][ jcoor ];  // \div d at this quadrature point ( - trace( A ) )

        for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
            A[ jcoor ][ jcoor ] += s;  // I\div d - (\grad d)^T at this quadrature point

        s = 0;
        for ( icoor = 0; icoor < fe_u.nbCoor(); icoor++ )
            for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
                s += guk[ icoor ][ jcoor ] * A[ icoor ][ jcoor ]; // \grad u^k : [I\div d - (\grad d)^T] at each quadrature point
        aux[ ig ] = s;
    }

    //
    // Numerical integration
    //

    // Loop on nodes, i.e. loop on elementary vector components
    for ( i = 0; i < fe_p.nbFEDof(); i++ )
    {

        // loop on quadrature points
        s = 0;
        for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
        {
            //std::cout << i << " " << ig << " " << fe_p.phi(i, ig) << " " << aux[ig] << std::endl;
            s += aux[ ig ] * fe_p.phi( i, ig ) * fe_u.weightDet( ig );
        }
        vec [ i ] += coef * s;
    }
}


//
// coef * ( [I\div d - (\grad d)^T - \grad d] \grap p, \grad q  ) for Newton FSI
//
void source_press2( Real coef, const ElemVec& p_loc, const ElemVec& d_loc, ElemVec& elvec,
                    const CurrentFE& fe, int iblock )
{

    Real A[ fe.nbCoor() ][ fe.nbCoor() ];     //  I\div d - (\grad d)^T - \grad d at a quadrature point
    Real B[ fe.nbCoor() ][ fe.nbCoor() ];     // - \grad d
    Real gpk[ fe.nbCoor() ];   // \grad p^k at a quadrature point
    Real aux[ fe.nbQuadPt() ][fe.nbCoor()];              // [I\div d - (\grad d)^T - \grad d ]\grad p^k at each quadrature point

    ElemVec::vector_view vec = elvec.block( iblock );

    Real s, sA, sG;
    UInt icoor, jcoor, ig, i;


    // loop on quadrature points
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {


        // loop on space coordinates
        for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            sG = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
                sG += fe.phiDer( i, icoor, ig ) * p_loc.vec() [ i ]; //  \grad p^k at a quadrature point
            gpk[ icoor ] = sG;

            // loop  on space coordinates
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                sA = 0.0;
                for ( i = 0; i < fe.nbFEDof(); i++ )
                    sA -= fe.phiDer( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe.nbFEDof() ]; //  - (\grad d) ^T at a quadrature point
                A[ icoor ][ jcoor ] = sA;
                B[ jcoor ][ icoor ] = sA;
            }
        }

        s = 0.0;
        for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            s -= A[ jcoor ][ jcoor ];  // \div d at this quadrature point ( - trace( A ) )

        for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            A[ jcoor ][ jcoor ] += s;  // I\div d - (\grad d)^T at this quadrature point

        for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                A[ icoor ][ jcoor ] += B[ icoor ][ jcoor ]; // I\div d - (\grad d)^T - \grad d at this quadrature point

        s = 0;
        for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                s += A[ icoor ][ jcoor ] * gpk[jcoor]; // [I\div d - (\grad d)^T -\grad d] \grad p^k at each quadrature point
            aux[ ig ][icoor] = s;
        }
    }

    //
    // Numerical integration
    //

    // Loop on nodes, i.e. loop on elementary vector components
    for ( i = 0; i < fe.nbFEDof(); i++ )
    {

        // loop on quadrature points
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                s += aux[ ig ][icoor] * fe.phiDer( i, icoor, ig ) * fe.weightDet( ig );
        vec [ i ] += coef * s;
    }
}





//
//Cholesky decomposition
void choldc( KNM<Real> &a, KN<Real> &p )
{
    int i, j, k;
    Real sum;

    int n = a.N();
    for ( i = 0; i < n; i++ )
    {
        for ( j = i; j < n; j++ )
        {
            for ( sum = a( i, j ), k = i - 1; k >= 0; k-- )
                sum -= a( i, k ) * a( j, k );
            if ( i == j )
            {
                p( i ) = sqrt( sum );
            }
            else
                a( j, i ) = sum / p( i );
        }
    }
}
//
//Cholesky solution
void cholsl( KNM<Real> &a, KN<Real> &p, KN<Real> &b, KN<Real> &x )
{
    int i, k;
    Real sum;

    int n = a.N();
    for ( i = 0; i < n; i++ )
    {
        for ( sum = b( i ), k = i - 1; k >= 0; k-- )
            sum -= a( i, k ) * x( k );
        x( i ) = sum / p( i );
    }
    for ( i = n - 1; i >= 0; i-- )
    {
        for ( sum = x( i ), k = i + 1; k < n; k++ )
            sum -= a( k, i ) * x( k );
        x( i ) = sum / p( i );
    }
}

//
// coef * (  (\grad u^k):[I\div d - (\grad d)^T] , q  ) for Newton FSI
//
void source_press( Real coef, const ElemVec& uk_loc, ElemMat& elmat,
                   const CurrentFE& fe_u, const CurrentFE& fe_p, ID mmDim , int iblock )
{

    Real A[ fe_u.nbCoor() ][ fe_u.nbCoor() ][fe_u.nbFEDof()][ fe_u.nbCoor() ];     //  I\div d - (\grad d)^T at a quadrature point
    Real guk[ fe_u.nbCoor() ][ fe_u.nbCoor() ];   // \grad u^k at a quadrature point
    Real aux[ fe_u.nbQuadPt() ][fe_u.nbFEDof()][ fe_u.nbCoor() ];              // grad u^k:[I\div d - (\grad d)^T] at each quadrature point

    //    double A[ fe_u.nbCoor() ][ fe_u.nbCoor() ][fe_u.nbFEDof()][ fe_u.nbCoor() ];

    Real /*s,*/ l/*[fe_u.nbFEDof()][fe_u.nbCoor()]*/, sA/*[fe_u.nbFEDof()][fe_u.nbCoor()]*/, sG;
    UInt icoor, jcoor, ig;


    // loop on quadrature points
    for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
    {


        ////INIT
        for (UInt p=0; p<fe_u.nbCoor(); ++p)
        {
            for (UInt q=0; q<fe_u.nbCoor(); ++q)
            {
                for (UInt d=0; d<fe_u.nbFEDof(); ++d)
                {
                    for (UInt e=0; e<fe_u.nbCoor(); ++e)
                        A[p][q][d][e]=0.;
                }
                //guk[p][q]=0.;
            }
            //for(int h=0; h<fe_u.nbFEDof(); ++h)
            //aux[ig][h][p]=0.;
        }
        ////END INIT

        // loop on space coordinates
        for ( icoor = 0; icoor < fe_u.nbCoor(); icoor++ )
        {

            //for ( i = 0;i < fe_u.nbFEDof();i++ )
            //for ( short k = 0;k < fe_u.nbCoor();k++ )
            //sA[i][k] = 0.0;
            // loop  on space coordinates
            for ( UInt kcoor = 0; kcoor < fe_u.nbCoor(); kcoor++ )
                for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
                {
                    sG = 0.0;
                    for ( UInt i = 0; i < fe_u.nbFEDof(); i++ )
                    {
                        sG += fe_u.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbFEDof() ]; //  \grad u^k at a quadrature point
                        sA = -fe_u.phiDer( i, icoor, ig )/**fe_u.phi( i, ig )*/ ;/** d_loc.vec() [ i + jcoor * fe_u.nbFEDof() ];*/ //  - (\grad d) ^T at a quadrature point
                        A[ icoor ][ jcoor ][i][jcoor] = sA; // -\delta_{jcoor kcoor} \partial_{icoor}
                    }
                    guk[ icoor ][ jcoor ] = sG;
                }
        }





        double z[fe_u.nbCoor()];
        for ( UInt i = 0; i < fe_u.nbFEDof(); i++ )
        {
            //                for ( int kcoor = 0;kcoor < fe_u.nbCoor();kcoor++ )
            for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
            {
                //if(icoor==jcoor)
                z[ jcoor ] = A[ jcoor ][ jcoor ][i][jcoor];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
            }

            for ( jcoor = 0; jcoor < fe_p.nbCoor(); jcoor++ )
            {
                for (UInt kcoor = 0; kcoor < fe_p.nbCoor(); kcoor++ )
                {
                    //if(icoor==jcoor)
                    A[ jcoor ][ jcoor ][i][kcoor] -= z[ kcoor ];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
                }
            }



//                         for ( jcoor = 0;jcoor < fe_u.nbCoor();jcoor++ )
//                             {
//                             //if(kcoor==jcoor)
//                                     for ( icoor = 0;icoor < fe_u.nbCoor();icoor++ )
//                                         {
//                                             //    if(icoor==jcoor)
//                                             A[ icoor ][ jcoor ][i][kcoor] -= A[ kcoor ][ icoor ][i][jcoor];  // \delta_{jcoor icoor}\partial_{kcoor}
//                                             //else
//                                             //f[icoor][jcoor]=0.;
//                                         }
//                             }
            //                    }
            //                for ( UInt kcoor = 0;kcoor < fe_u.nbCoor();kcoor++ )
            //                    {
//                         for ( jcoor = 0;jcoor < fe_u.nbCoor();jcoor++ )
//                         {
//                             for ( icoor = 0;icoor < fe_u.nbCoor();icoor++ )
//                                 {
//                                     //for ( jcoor = 0;jcoor < fe_u.nbCoor();jcoor++ )
//                                     if(icoor==jcoor)
//                                         A[ icoor ][ jcoor ][i][kcoor] += f[icoor][jcoor];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
//                 //                for ( i = 0;i < fe_u.nbFEDof();i++ )
//                 //                    for(jcoor=0;jcoor<fe_u.nbCoor();++jcoor)
//                 //                        s[i][jcoor] = 0.0;
//                               //for ( jcoor = 0;jcoor < fe_u.nbCoor();jcoor++ )
//                               //{
//                                 }
//                         }
            //                    }
            //                for ( UInt kcoor = 0;kcoor < fe_u.nbCoor();kcoor++ )
            //                    {
            for ( UInt kcoor = 0; kcoor < fe_u.nbCoor(); kcoor++ )
            {
                l=0;
                //for ( kcoor = 0;kcoor < fe_u.nbCoor();kcoor++ )
                for ( icoor = 0; icoor < fe_u.nbCoor(); icoor++ )
                    for ( jcoor = 0; jcoor < fe_u.nbCoor(); jcoor++ )
                    {
                        l += guk[ icoor ][ jcoor ] * A[ icoor ][ jcoor ][i][kcoor]; // \grad u^k : [I\div d - (\grad d)^T] at each quadrature point
                        //}
                    }
                aux[ ig ][i][kcoor] = l;
                //}
            }
            //std::cout<<"aux[ ig ][i][jcoor] "<<  aux[ ig ][i][jcoor] <<std::endl;
        }
    }
    //
    // Numerical integration
    //
    for ( UInt kcoor = 0; kcoor < fe_u.nbCoor(); kcoor++ )
    {
        //            for ( short l = 0;i < fe_u.nbFEDof();i++ )
        //for ( jcoor = 0;jcoor < fe_u.nbCoor();jcoor++ )
        double l = 0.;

        ElemMat::matrix_view mat = elmat.block( iblock, kcoor );
        for ( UInt j = 0; j < mmDim; j++ )
            for ( UInt i = 0; i < fe_p.nbFEDof(); i++ )
                mat(i,j)=0.;

        // Loop on nodes, i.e. loop on elementary vector components
        for ( UInt j = 0; j < mmDim; j++ )
            for (UInt i = 0; i < fe_p.nbFEDof(); i++ )
            {
                l=0.;
                // loop on quadrature points
                for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
                {
                    //            std::cout << ig << " " << fe_p.phi(i, ig) << std::endl;
                    l += aux[ ig ][j][kcoor] * fe_p.phi( i, ig ) * fe_u.weightDet( ig );
                }
                mat( i , j) += coef * l;
                //std::cout<<"mat [ i , j] "<<mat [ i , j]<<std::endl;
                //std::cout<<"l "<<l<<std::endl;
            }
    }
}



//shape_terms_vel:
//source_mass1
// coef * ( - \grad w^k :[I\div d - (\grad d)^T] u^k + convect^T[I\div d - (\grad d)^T] (\grad u^k)^T , v  ) for Newton FSI
//
// Remark: convect = u^n-w^k relative vel.
//
//source_stress
// coef * ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  ) for Newton FSI
//
//source_stress2
// \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )
//
//source_mass3
// 0.5 * ( \grad u^n :[2 I \div d - (\grad d)^T]  u^k , v  ) for Newton FSI
//
//optionally (if full implicit):
//source_mass2
// -rho * ( \grad u^k dw, v  )
//
//convective term
// rho * ( \grad u^k du, v  )
//
void shape_terms(
    //const ElemVec& d_loc,
    Real rho,
    Real mu,
    const ElemVec& un_loc,
    const ElemVec& uk_loc,
    const ElemVec& wk_loc,
    const ElemVec& convect_loc,
    const ElemVec& pk_loc,
    ElemMat& elmat,
    const CurrentFE& fe,
    const CurrentFE& fe_p,
    ID /*mmDim*/,
    ElemMat& /*elmatP*/,
    int /*iblock*/,
    bool wImplicit,
    Real alpha,
    boost::shared_ptr<ElemMat> elmat_convect
)
{
    //Real BGrConv[ fe.nbCoor() ][ fe.nbCoor() ]/*[fe.nbFEDof()]*/;                 // \grad (convect) at a quadrature point
    Real eta[ fe.nbCoor() ][ fe.nbCoor() ][fe.nbFEDof()][ fe.nbCoor() ];                 // I\div d - (\grad d)^T at a quadrature point
    Real etaMass3[ fe.nbCoor() ][ fe.nbCoor() ][fe.nbFEDof()][ fe.nbCoor() ];                 // I\div d - (\grad d)^T at a quadrature point
    Real uk[ fe.nbQuadPt() ][ fe.nbCoor() ];              // u^k quadrature points
    Real guk[ fe.nbQuadPt() ][ fe.nbCoor() ][ fe.nbCoor() ];  // \grad u^k at quadrature points
    Real duk[ fe.nbQuadPt() ];  // \div u^k at quadrature points
    Real gun[ fe.nbQuadPt() ][ fe.nbCoor() ][ fe.nbCoor() ];  // \grad u^k at quadrature points
    Real convect[ fe.nbCoor() ];                      // convect at quadrature points
    Real convect_eta[ fe.nbQuadPt() ][ fe.nbCoor() ][fe.nbFEDof()][ fe.nbCoor() ];       // (convect)^T [I\div d - (\grad d)^T] at quadrature points
    Real sigma[ fe.nbCoor() ][ fe.nbCoor() ];             // [-p^k I + 2*mu e(u^k)] a quadrature point
    Real B[ fe.nbQuadPt() ][ fe.nbCoor() ][ fe.nbCoor() ][fe.nbFEDof()][ fe.nbCoor() ];  // [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] at each quadrature point
    Real gd[ fe.nbCoor() ][ fe.nbCoor() ][fe.nbFEDof()][ fe.nbCoor() ];                // \grad d at a quadrature point
    Real A2[ fe.nbQuadPt() ][ fe.nbCoor() ][ fe.nbCoor() ][fe.nbFEDof()][ fe.nbCoor() ];  // \grad u^k \grad d + [\grad d]^T[\grad u^k]^T  at each quadrature point
    Real aux[ fe.nbQuadPt() ][ fe.nbCoor() ][fe.nbFEDof()][fe.nbCoor()];
    Real BMass[ fe.nbCoor() ][ fe.nbCoor() ];
    Real auxMass[ fe.nbQuadPt() ][fe.nbFEDof()][fe.nbCoor()];
    Real auxMass3[ fe.nbQuadPt() ][fe.nbFEDof()][fe.nbCoor()];

    Real s,  sA, sB, sGk, sGn, pk;

    UInt icoor, jcoor, ig, kcoor, i, j;

    // loop on quadrature points

    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        ////INIT
        for (UInt p=0; p<fe.nbCoor(); ++p)
        {
            uk[ig][p]=0.;
            convect[p]=0.;
            for (UInt q=0; q<fe.nbCoor(); ++q)
            {
                guk[ig][p][q]=0.;
                gun[ig][p][q]=0.;
                sigma[p][q]=0.;
                BMass[p][q]=0.;
                for (UInt d=0; d<fe.nbFEDof(); ++d)
                {
                    auxMass[ ig ][d][q]=0.;
                    auxMass3[ ig ][d][q]=0.;
                    aux[ig][p][d][q]=0.;
                    convect_eta[ig][p][d][q]=0.;
                    for (UInt e=0; e<fe.nbCoor(); ++e)
                    {
                        gd[p][q][d][e]=0.;
                        eta[p][q][d][e]=0.;
                        etaMass3[p][q][d][e]=0.;
                        A2[ig][p][q][d][e]=0.;
                        B[ig][p][q][d][e]=0.;
                    }
                }
            }
        }
        ////////END INIT

        // loop on space coordindates
        for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {

            // each compontent of uk at each quadrature points
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
                s += fe.phi( i, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
            uk[ ig ][ icoor ] = s;//uk_x(pt_ig), uk_y(pt_ig), uk_z(pt_ig)

            // each compontent of convect at this quadrature point
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
                s += fe.phi( i, ig ) * convect_loc.vec() [ i + icoor * fe.nbFEDof() ];
            convect[ icoor ] = s;


            // loop  on space coordindates
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                sB = 0.0;
                sGk = 0.0;
                sGn = 0.0;
                for ( i = 0; i < fe.nbFEDof(); i++ )
                {
                    gd[ icoor ][ jcoor ][i][jcoor] = fe.phiDer( i, jcoor, ig )/** d_loc.vec() [ i + jcoor * fe.nbFEDof() ]*/;

                    sGk += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at each quadrature point
                    sGn += fe.phiDer( i, jcoor, ig ) * un_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at each quadrature point
                    sB -= fe.phiDer( i, jcoor, ig ) * wk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad (- w^k) at this quadrature point
                    sA = -fe.phiDer( i, icoor, ig ) /** d_loc.vec() [ i + jcoor * fe.nbFEDof() ]*/; //  - (\grad d) ^T at this quadrature point
                    eta[ icoor ][ jcoor ][i][ jcoor ] = sA; // -(\grad d) ^T at this quadrature point
                    etaMass3[ icoor ][ jcoor ][i][ jcoor ] = sA; // -(\grad d) ^T at this quadrature point
                }
                guk[ ig ][ icoor ][ jcoor ] = sGk; // \grad u^k at each quadrature point
                gun[ ig ][ icoor ][ jcoor ] = sGn; // \grad u^n at each quadrature point
                BMass[ icoor ][ jcoor ] = sB;
            }
        }

        //!a part of source_stress
        pk = 0.0;
        for ( i = 0; i < fe_p.nbFEDof(); i++ )
            pk += fe_p.phi( i, ig ) * pk_loc.vec() [ i ]; // p^k at this quadrature point

        // sigma = [-p^k I + 2*mu e(u^k)] a quadrature point
        for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
        {
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                sigma[ icoor ][ jcoor ] = mu * ( guk[ig][ icoor ][ jcoor ] + guk[ig][ jcoor ][ icoor ] );
            }
            sigma[ icoor ][ icoor ] -= pk;
        }

        //!building the tensor \f$\eta = [I\nabla\cdot d - (\nabla d)^T]\f$
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {

            Real z[fe.nbCoor()];
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                //if(icoor==jcoor)
                z[ jcoor ] = eta[ jcoor ][ jcoor ][i][jcoor];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
            }
            for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
            {
                for ( kcoor = 0; kcoor < fe.nbCoor(); kcoor++ )
                {
                    eta[ jcoor ][ jcoor ][i][kcoor] -= z[kcoor];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
                }
            }

            //!source_mass1

            for ( kcoor = 0; kcoor < fe.nbCoor(); kcoor++ )
            {
                s = 0;
                for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                    for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                        s += BMass[ icoor ][ jcoor ] * eta[ icoor ][ jcoor ][i][kcoor]; // \grad (-w^k):[I\div d - (\grad d)^T] at each quadrature point
                auxMass[ ig ][i][kcoor] = s;

                for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                {
                    s = 0.;
                    for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                    {
                        s += convect[ icoor ] * eta[ icoor ][ jcoor ][i][kcoor]; // convect^T [I\div d - (\grad d)^T]
                    }
                    convect_eta[ ig ][ jcoor ][i][kcoor] = s;
                }
            }
            //! source_stress

            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
            {
                for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                {
                    for ( kcoor = 0; kcoor < fe.nbCoor(); kcoor++ )
                    {
                        s = 0.;
                        for (UInt zcoor = 0; zcoor < fe.nbCoor(); zcoor++ )
                            s += sigma[ icoor ][ zcoor ] * eta[ zcoor ][ jcoor ][i][kcoor];
                        B[ ig ][ icoor ][ jcoor ][i][kcoor] = s;
                    }
                }
            }


            //! source_stress2
            for ( kcoor = 0; kcoor < fe.nbCoor(); kcoor++ )
            {
                for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                {
                    for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                    {
                        s = 0.;
                        for ( UInt zcoor = 0; zcoor < fe.nbCoor(); zcoor++ )
                            // \grad u^k \grad d + [\grad d]^T[\grad u^k]^T  at each quadrature point
                            s += guk[ig][ icoor ][ zcoor ] * gd[ zcoor ][ jcoor ][i][kcoor] + gd[ zcoor ][ icoor ][i][kcoor] * guk[ig][ jcoor ][ zcoor ];
                        A2[ ig ][ icoor ][ jcoor ][i][kcoor] = s;
                    }
                }
            }

            if (wImplicit)//source_mass2 and convective term derivative
            {
                for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                {
                    //s = 0.0;
                    for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                    {
                        aux[ ig ][ icoor ][i][jcoor] = guk[ig][ icoor ][ jcoor ]  * fe.phi( i, ig ) /**alpha*/;
                    }

                }
                for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                {
                    duk[ig] += guk[ig][ jcoor ][ jcoor ];  //div uk
                }
            }


            //source_mass3
            // coef * ( \grad u^n :[2 I \div d - (\grad d)^T]  u^k , v  ) for Newton FSI
            //
            //

            //!building the tensor \f$\eta = [I\nabla\cdot d - (\nabla d)^T]\f$
            //                for ( int kcoor = 0;kcoor < fe.nbCoor();kcoor++ )

            for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                {
                    for (kcoor = 0; kcoor < fe.nbCoor(); kcoor++ )
                    {
                        //if(icoor==jcoor)
                        etaMass3[ icoor ][ jcoor ][i][kcoor] -= 2*z[kcoor];  //! -2*\delta_{jcoor, kcoor} \partial_{icoor} + \delta_{jcoor, icoor}\partial_{kcoor}
                    }
                }

            for (kcoor =0; kcoor<fe.nbCoor(); ++kcoor)
            {
                s = 0;
                for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
                    for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                        s += gun[ig][ icoor ][ jcoor ] * etaMass3/**/[ icoor ][ jcoor ][i][kcoor]; // \grad u^n:[/*2*/ * I\div d - (\grad d)^T] at each quadrature point
                auxMass3[ ig ][i][kcoor] = s;
                //if fullImplicit un=uk
            }

            ///////////////////////////////////////////////////////////////////////
        }

    }



    //
    // Numerical integration
    //
    Real g=0.;
    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0; icoor < fe.nbCoor(); icoor++ )
    {
        for ( kcoor = 0; kcoor < fe.nbCoor(); kcoor++ )
        {

            // the block iccor of the elementary vector
            ElemMat::matrix_view mat = elmat.block( icoor, kcoor );
            boost::shared_ptr<ElemMat::matrix_view> mat_convect;

            if (elmat_convect.get())
            {
                mat_convect.reset(new ElemMat::matrix_view(elmat_convect->block( icoor, kcoor )));
            }
            // loop on nodes, i.e. loop on components of this block
            for ( i = 0; i < fe.nbFEDof(); i++ )
            {
                for ( j = 0; j < fe.nbFEDof(); j++ )
                {

                    // loop on quadrature points
                    s = 0.;
                    g = 0.;

                    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
                    {
                        // - convect^T [I\div d - (\grad d)^T] (u^k)^T :\grad\phi_i
                        for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                        {
                            //source_mass1
                            s += convect_eta[ ig ][ jcoor ][j][kcoor] * guk[ ig ][ icoor ][jcoor] * fe.phi( i, ig ) * fe.weightDet( ig )*rho;
                            //source_stress
                            s += B[ ig ][ icoor ][ jcoor ][j][kcoor] * fe.phiDer( i, jcoor, ig ) * fe.weightDet( ig );
                            //source_stress2
                            s -= fe.phiDer( i, jcoor, ig ) * A2[ ig ][ icoor ][ jcoor ][j][kcoor] * fe.weightDet( ig )*mu;


                        }
                        //source_mass1
                        s += auxMass[ ig ][j][kcoor] * uk[ ig ][ icoor ] * fe.phi( i, ig ) * fe.weightDet( ig )*rho;

                        //source_mass3
                        s += 0.5*auxMass3[ ig ][j][kcoor] * uk[ ig ][ icoor ] * fe.phi( i, ig ) * fe.weightDet( ig )*rho;

                        if (wImplicit)
                        {
                            //source_mass2{
                            s -= aux[ ig ][ icoor ][j][kcoor] * fe.phi( i, ig ) * fe.weightDet( ig )*rho*alpha;
                            //convective term
                            g += aux[ ig ][ icoor ][j][kcoor] * fe.phi( i, ig ) * fe.weightDet( ig )*rho;
                        }
                    }
                    mat( i , j ) += s;
                    if (wImplicit && mat_convect.get())
                        (*mat_convect)( i , j ) += g;
                }
            }
        }
    }
}


//----------------------------------------------------------------------
// Compute the gradient in the Hdiv space, i.e. the opposite and transpose of the divergence matrix.
void grad_Hdiv( Real coef, ElemMat& elmat, const CurrentFE& dualFE, const CurrentFE& primalFE, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Real sumdivphi(0.);
    const QuadRule& dualQuadRule( dualFE.quadRule() );

    // Loop over all the degrees of freedom of the dual variable.
    for ( UInt i(0); i < dualFE.nbFEDof(); ++i )
    {
        // Loop over all the degrees of freedom of the primal variable.
        for ( UInt j(0); j < primalFE.nbFEDof(); ++j )
        {
            sumdivphi = 0.;
            // Loop over all the quadrature points.
            for ( UInt ig(0); ig < dualFE.nbQuadPt(); ++ig )
            {
                // There is no jacobian because of property of Piola transformation.
                sumdivphi -= primalFE.phi( j, ig ) * dualFE.divPhiRef( i, ig ) * dualQuadRule.weight( ig );
            }
            mat( i, j ) += coef * sumdivphi;
        }
    }
}

//----------------------------------------------------------------------
// Compute the divergence in the Hdiv space.
void div_Hdiv( Real coef, ElemMat& elmat, const CurrentFE& dualFE, const CurrentFE& primalFE, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Real sumdivphi(0.);
    const QuadRule& dualQuadRule( dualFE.quadRule() );

    // Loop over all the degrees of freedom of the dual variable.
    for ( UInt i(0); i < dualFE.nbFEDof(); ++i )
    {
        // Loop over all the degrees of freedom of the primal variable.
        for ( UInt j(0); j < primalFE.nbFEDof(); ++j )
        {
            sumdivphi = 0.;
            // Loop over all the quadrature points.
            for ( UInt ig(0); ig < dualFE.nbQuadPt(); ++ig )
            {
                // There is no jacobian because of property of Piola transformation.
                sumdivphi += primalFE.phi( j, ig ) * dualFE.divPhiRef( i, ig ) * dualQuadRule.weight( ig );
            }
            mat( j, i ) += coef * sumdivphi;
        }
    }
}

//----------------------------------------------------------------------
// Compute a Hdiv function dot product with the outwart unit normal times a hybrid function.
void TP_VdotN_Hdiv( Real coef, ElemMat& elmat, const RefFEHybrid& hybridFE,
                    const RefFEHybrid& dualDotNFE, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt nbnode;
    Real sum(0.);

    // Loop over all the staticBdFE.
    for ( UInt nf(0); nf < hybridFE.numberBoundaryFE(); ++nf )
    {
        // Take the staticBdFE of the hybrid finite element.
        const StaticBdFE & boundaryElementHybridFE( hybridFE[ nf ] );
        // Take the staticBdFE of the Hdiv function dot product outward unit normal.
        const StaticBdFE & boundaryElementDualDotNFE( dualDotNFE[ nf ] );
        nbnode = boundaryElementHybridFE.nbNode();

        // Loop over all the the degrees of freedom of the dual dot normal variable.
        for ( UInt i(0); i < nbnode; ++i )
        {
            // Loop over all the degrees of freedom of the hybrid variable.
            for ( UInt j(0); j < nbnode; ++j )
            {
                sum = 0.;
                // Loop over all the quadrature point.
                for ( UInt ig(0); ig < boundaryElementHybridFE.nbQuadPt(); ++ig )
                    // Using the Piola transform properties.
                    sum += boundaryElementHybridFE.phi( j , ig ) *
                           boundaryElementDualDotNFE.phi( i , ig ) *
                           boundaryElementHybridFE.weightMeas( ig );

                // The matrix is block diagonal, so the size of the blocks is bdfe.nbNode.
                mat( nf * nbnode + i, nf * nbnode + j ) += sum * coef;
            }
        }
    }
}

//----------------------------------------------------------------------
// Compute the mass matrix for hybrid variable.
void TP_TP_Hdiv( Real coef, ElemMat& elmat, const RefFEHybrid& hybridFE, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt nbnode;
    Real sum(0.);

    // Loop over all the staticBdFE.
    for ( UInt nf(0); nf < hybridFE.numberBoundaryFE(); ++nf )
    {
        // Take the staticBdFE of the hybrid finite element.
        const StaticBdFE & boundaryElementHybridFE( hybridFE[ nf ] );
        nbnode = boundaryElementHybridFE.nbNode();

        // Loop over all the degrees of freedom of the first hybrid variable.
        for ( UInt i(0); i < nbnode; ++i )
        {
            // Loop over all the the degrees of freedom of the second hybrid variable.
            for ( UInt j(0); j < nbnode; ++j )
            {
                sum = 0.;
                // Loop over all the quadrature point.
                for ( UInt ig(0); ig < boundaryElementHybridFE.nbQuadPt() ; ++ig )
                    // Using the Piola transform properties.
                    sum += boundaryElementHybridFE.phi( j , ig ) *
                           boundaryElementHybridFE.phi( i , ig ) *
                           boundaryElementHybridFE.weightMeas( ig );

                // The matrix is block diagonal, so the size of the blocks is bdfe.nbNode.
                mat( nf * nbnode + i, nf * nbnode + j ) += sum * coef;
            }
        }
    }
}


//----------------------------------------------------------------------
// Compute the mass matrix in Hdiv with a real scalar coefficient.
void mass_Hdiv( Real coef, ElemMat& elmat, const CurrentFE& dualFE, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Real sum(0.);

    // Loop over all the degrees of freedom of the first dual variable.
    for ( UInt j(0); j < dualFE.nbFEDof(); ++j )
    {
        // Loop over all the degrees of freedom of the second dual variable.
        for ( UInt i(0); i < dualFE.nbFEDof() /* by symmetry j+1 */; ++i )
        {
            sum = 0.;
            // Loop over all the quadrature point.
            for ( UInt ig(0); ig < dualFE.nbQuadPt(); ++ig )
            {
                // Loop over all the space dimension, e.g. in 3D three times, in 2D two times.
                for ( UInt icoor(0); icoor < dualFE.nbCoor() ; ++icoor )
                {
                    sum += dualFE.phi( j, icoor, ig ) * dualFE.phi( i, icoor, ig ) * dualFE.wDetJacobian( ig );
                }
            }
            // Beware coef is the inverse of the permeability, not the permeability.
            mat( i, j ) += sum * coef;
        }
    }
}

//----------------------------------------------------------------------
// Compute the mass matrix in Hdiv with a real matrix coefficient.
void mass_Hdiv( Matrix const&  Invperm, ElemMat& elmat, const CurrentFE& dualFE, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Real partialSum(0.), sum(0.);

    // Loop over all the degrees of freedom of the first dual variable.
    for ( UInt j(0); j < dualFE.nbFEDof(); ++j )
    {
        // Loop over all the degrees of freedom of the second dual variable.
        for ( UInt i(0); i < dualFE.nbFEDof() /* by symmetry j+1 */; ++i )
        {
            sum = 0.;
            // Loop over all the quadrature point.
            for ( UInt ig(0); ig < dualFE.nbQuadPt(); ++ig )
            {
                // partialSum = phi[icoor]^T * K^-1 * phi[jcoor] for all icorr and jcoor.
                partialSum = 0.;
                /* Loop over all the space dimension of the first dual variable,
                    e.g. in 3D three times, in 2D two times.*/
                for ( UInt icoor(0); icoor < dualFE.nbCoor(); ++icoor )
                {
                    /* Loop over all the space dimension of the first dual variable,
                          e.g. in 3D three times, in 2D two times.*/
                    for ( UInt jcoor(0); jcoor < dualFE.nbCoor(); ++jcoor )
                    {
                        // Invperm is the inverse of the permeability.
                        partialSum += ( Invperm( icoor, jcoor ) *
                                        dualFE.phi( j, jcoor, ig ) *
                                        dualFE.phi( i, icoor, ig ) );
                    }
                }
                sum += partialSum * dualFE.wDetJacobian( ig );
            }
            mat( i, j ) += sum ;
        }
    }
}

//----------------------------------------------------------------------
// Compute the mass matrix in Hdiv with a real function coefficient.
void mass_Hdiv( Real ( *InvpermFun ) ( const Real&, const Real&, const Real& ),
                ElemMat& elmat, const CurrentFE& dualFE, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Real sum(0.), x(0.), y(0.), z(0.);

    // Loop over all the degrees of freedom of the first dual variable.
    for ( UInt j(0); j < dualFE.nbFEDof(); ++j )
    {
        // Loop over all the degrees of freedom of the second dual variable.
        for ( UInt i(0); i < dualFE.nbFEDof() /* by symmetry j+1 */; ++i )
        {
            sum = 0.;
            // Loop over all the quadrature point.
            for ( UInt ig(0); ig < dualFE.nbQuadPt(); ++ig )
            {
                // Get the coordinate in the current element of the current quadrature point.
                x = dualFE.quadNode(ig, 0);
                y = dualFE.quadNode(ig, 1);
                z = dualFE.quadNode(ig, 2);

                // Loop over all the space dimension, e.g. in 3D three times, in 2D two times.
                for ( UInt icoor(0); icoor < dualFE.nbCoor(); ++icoor )
                {
                    // Caution Invperm is the inverse of the permeability.
                    sum += InvpermFun( x, y, z ) *
                           dualFE.phi( j, icoor, ig ) *
                           dualFE.phi( i, icoor, ig ) *
                           dualFE.wDetJacobian( ig );
                }
            }
            mat( i, j ) += sum;
        }
    }
}



}

#endif
