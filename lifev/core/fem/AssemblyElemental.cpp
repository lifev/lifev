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

#include <lifev/core/fem/AssemblyElemental.hpp>
#include <boost/multi_array.hpp>

namespace LifeV
{

namespace AssemblyElemental
{
void mass (MatrixElemental& localMass,
           const CurrentFE& massCFE,
           const Real& coefficient,
           const UInt& fieldDim)
{
    const UInt nbFEDof (massCFE.nbFEDof() );
    const UInt nbQuadPt (massCFE.nbQuadPt() );
    Real localValue (0);

    MatrixElemental::matrix_type mat_tmp (nbFEDof, nbFEDof);

    // Loop over the basis functions (diagonal part)
    for (UInt iDof (0); iDof < nbFEDof; ++iDof)
    {
        localValue = 0.0;
        for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
        {
            localValue += massCFE.phi (iDof, iQuadPt)
                          * massCFE.phi (iDof, iQuadPt)
                          * massCFE.wDetJacobian (iQuadPt);
        }
        localValue *= coefficient;

        // Add on the local matrix
        mat_tmp (iDof, iDof) = localValue;
    }

    // Loop over the basis functions (extradiagonal part)
    for (UInt iDof (0); iDof < nbFEDof; ++iDof)
    {
        // Build the local matrix only where needed:
        // Lower triangular + diagonal parts
        for (UInt jDof (0); jDof < iDof; ++jDof)
        {
            localValue = 0.0;

            //Loop on the quadrature nodes
            for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
            {
                localValue += massCFE.phi (iDof, iQuadPt)
                              * massCFE.phi (jDof, iQuadPt)
                              * massCFE.wDetJacobian (iQuadPt);
            }

            localValue *= coefficient;

            // Add on the local matrix
            mat_tmp (iDof, jDof) = localValue;
            mat_tmp (jDof, iDof) = localValue;
        }
    }

    // Copying the mass in all the diagonal blocks (just one for scalar problem)
    for (UInt iDim (0); iDim < fieldDim; ++iDim)
    {
        MatrixElemental::matrix_view mat = localMass.block (iDim, iDim);
        mat += mat_tmp;
    }
}

void stiffness (MatrixElemental& localStiff,
                const CurrentFE& stiffCFE,
                const Real& coefficient,
                const UInt& fieldDim)
{
    const UInt nbFEDof (stiffCFE.nbFEDof() );
    const UInt nbQuadPt (stiffCFE.nbQuadPt() );
    Real localValue (0);

    MatrixElemental::matrix_type mat_tmp (nbFEDof, nbFEDof);

    // Loop over the basis functions (diagonal part)
    for (UInt iDof (0); iDof < nbFEDof; ++iDof)
    {
        localValue = 0.0;

        // Loop on the quadrature noodes
        for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
        {
            for (UInt iDim (0); iDim < stiffCFE.nbLocalCoor(); ++iDim)
            {
                localValue += stiffCFE.dphi (iDof, iDim, iQuadPt)
                              * stiffCFE.dphi (iDof, iDim, iQuadPt)
                              * stiffCFE.wDetJacobian (iQuadPt);
            }
        }
        localValue *= coefficient;

        // Add on the local matrix
        mat_tmp (iDof, iDof) = localValue;
    }

    // Loop over the basis functions (extradiagonal part)
    for (UInt iDof (0); iDof < nbFEDof ; ++iDof)
    {
        for (UInt jDof (0); jDof < iDof; ++jDof)
        {
            localValue = 0.0;

            // Loop on the quadrature nodes
            for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
            {
                for (UInt iDim (0); iDim < stiffCFE.nbLocalCoor(); ++iDim)
                {
                    localValue += stiffCFE.dphi (iDof, iDim, iQuadPt)
                                  * stiffCFE.dphi (jDof, iDim, iQuadPt)
                                  * stiffCFE.wDetJacobian (iQuadPt);
                }
            }


            localValue *= coefficient;

            // Put in the local matrix
            mat_tmp (iDof, jDof) = localValue;
            mat_tmp (jDof, iDof) = localValue;
        }
    }

    // Copying the diffusion in all the diagonal blocks
    for ( UInt iDim (0); iDim < fieldDim; ++iDim)
    {
        MatrixElemental::matrix_view mat = localStiff.block (iDim, iDim);
        mat += mat_tmp;
    }
}

void advectionNewton ( Real coef, VectorElemental& vel,
                       MatrixElemental& elmat, const CurrentFE& fe,
                       int iblock, int jblock )
{
    Real sum, sumDerivative;
    MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock, jblock );

    //Assemble the local matrix
    for ( UInt i = 0; i < fe.nbFEDof(); i++ )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); j++ )
        {
            sum = 0.;
            for ( UInt iq = 0; iq < fe.nbQuadPt(); iq++ )
            {
                //evaluate derivative
                VectorElemental::vector_view velicoor = vel.block ( iblock );
                sumDerivative = 0.;
                for ( UInt k = 0; k < fe.nbFEDof(); k++ )
                {
                    sumDerivative += velicoor ( k ) * fe.dphi ( k, jblock, iq );
                }
                //end evaluate derivative

                sum += fe.phi (j, iq) * sumDerivative * fe.phi (i, iq) * fe.weightDet ( iq );
            }
            mat_icomp ( i, j ) += sum * coef;
        }
    }
}

void grad ( MatrixElemental& elmat,
            const CurrentFE& uCFE,
            const CurrentFE& pCFE,
            const UInt& fieldDim)
{
    const UInt nbPFEDof (pCFE.nbFEDof() );
    const UInt nbUFEDof (uCFE.nbFEDof() );
    const UInt nbQuadPt (pCFE.nbQuadPt() );
    Real localValue (0);

    for (UInt iFDim (0); iFDim < fieldDim; ++iFDim)
    {
        MatrixElemental::matrix_view localView = elmat.block (iFDim, 0);

        for (UInt iDofP (0); iDofP < nbPFEDof; ++iDofP)
        {
            for (UInt iDofU (0); iDofU < nbUFEDof; ++iDofU)
            {
                localValue = 0.0;
                for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
                {
                    localValue += uCFE.dphi (iDofU, iFDim, iQuadPt)
                                  * pCFE.phi (iDofP, iQuadPt)
                                  * uCFE.wDetJacobian (iQuadPt);
                }
                localView (iDofU, iDofP) -= localValue;
            }
        }
    }
}

void divergence ( MatrixElemental& elmat,
                  const CurrentFE& uCFE,
                  const CurrentFE& pCFE,
                  const UInt& fieldDim,
                  const Real& coefficient)
{
    const UInt nbPFEDof (pCFE.nbFEDof() );
    const UInt nbUFEDof (uCFE.nbFEDof() );
    const UInt nbQuadPt (pCFE.nbQuadPt() );
    Real localValue (0);

    for (UInt iFDim (0); iFDim < fieldDim; ++iFDim)
    {
        MatrixElemental::matrix_view localView = elmat.block (0, iFDim);

        for (UInt iDofP (0); iDofP < nbPFEDof; ++iDofP)
        {
            for (UInt iDofU (0); iDofU < nbUFEDof; ++iDofU)
            {
                localValue = 0.0;
                for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
                {
                    localValue += uCFE.dphi (iDofU, iFDim, iQuadPt)
                                  * pCFE.phi (iDofP, iQuadPt)
                                  * uCFE.wDetJacobian (iQuadPt);
                }
                localView (iDofP, iDofU) -= coefficient * localValue;
            }
        }
    }
}

void stiffStrain (MatrixElemental& localStiff,
                  const CurrentFE& stiffCFE,
                  const Real& coefficient,
                  const UInt& fieldDim)
{
    const UInt nbFEDof (stiffCFE.nbFEDof() );
    const UInt nbQuadPt (stiffCFE.nbQuadPt() );
    Real localValue (0);
    const Real newCoefficient (coefficient * 0.5);

    stiffness (localStiff, stiffCFE, newCoefficient, fieldDim); // for the stiff part, we exploit the existing routine

    for ( UInt iFDim (0); iFDim < fieldDim; ++iFDim )
    {
        for ( UInt jFDim (0); jFDim < fieldDim; ++jFDim )
        {
            MatrixElemental::matrix_view localView = localStiff.block ( iFDim, jFDim );

            for ( UInt iDof (0); iDof < nbFEDof; ++iDof )
            {
                for ( UInt jDof (0); jDof < nbFEDof; ++jDof )
                {
                    localValue = 0.0;
                    for ( UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt )
                    {
                        localValue += stiffCFE.dphi ( iDof, jFDim, iQuadPt )
                                      * stiffCFE.dphi ( jDof, iFDim, iQuadPt )
                                      * stiffCFE.wDetJacobian ( iQuadPt );
                    }
                    localView ( iDof, jDof ) += newCoefficient * localValue;
                }
            }
        }
    }
}

void bodyForces (VectorElemental& localForce,
                 const CurrentFE& massRhsCFE,
                 const function_Type& fun,
                 const Real& t,
                 const UInt& fieldDim)
{

    const UInt nbFEDof (massRhsCFE.nbFEDof() );
    const UInt nbQuadPt (massRhsCFE.nbQuadPt() );
    std::vector<Real> fValues (nbQuadPt, 0.0);
    Real localValue (0.0);

    // Assemble the local diffusion
    for (UInt iterFDim (0); iterFDim < fieldDim; ++iterFDim)
    {
        VectorElemental::vector_view localView = localForce.block (iterFDim);

        // Compute the value of f in the quadrature nodes
        for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
        {
            fValues[iQuadPt] = fun (t,
                                    massRhsCFE.quadNode (iQuadPt, 0),
                                    massRhsCFE.quadNode (iQuadPt, 1),
                                    massRhsCFE.quadNode (iQuadPt, 2),
                                    iterFDim);
        }

        // Loop over the basis functions
        for (UInt iDof (0); iDof < nbFEDof ; ++iDof)
        {
            localValue = 0.0;

            //Loop on the quadrature nodes
            for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
            {
                localValue += fValues[iQuadPt]
                              * massRhsCFE.phi (iDof, iQuadPt)
                              * massRhsCFE.wDetJacobian (iQuadPt);
            }

            // Add on the local matrix
            localView (iDof) = localValue;
        }
    }

}

} // Namespace AssemblyElemental

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
void mass ( Real coef, MatrixElemental& elmat, const CurrentFE& fe,
            int iblock, int jblock )
/*
  Mass matrix: \int v_i v_j
*/
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    UInt i, ig;
    int iloc, jloc;
    Real s, coef_s;
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += fe.phi ( iloc, ig ) * fe.phi ( iloc, ig ) * fe.weightDet ( ig );
        }
        mat ( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += fe.phi ( iloc, ig ) * fe.phi ( jloc, ig ) * fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat ( iloc, jloc ) += coef_s;
        mat ( jloc, iloc ) += coef_s;
    }
}
//
// coeff*Mass
//
void mass ( Real coef, MatrixElemental& elmat, const CurrentFE& fe,
            int iblock, int jblock, UInt nb )
/*
  Mass matrix: \int v_i v_j (nb blocks on the diagonal, nb>1)
*/
{
    Matrix mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );
    UInt i, ig;
    int iloc, jloc;
    Real s, coef_s;
    mat_tmp = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += fe.phi ( iloc, ig ) * fe.phi ( iloc, ig ) * fe.weightDet ( ig );
        }
        mat_tmp ( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += fe.phi ( iloc, ig ) * fe.phi ( jloc, ig ) * fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat_tmp ( iloc, jloc ) += coef_s;
        mat_tmp ( jloc, iloc ) += coef_s;
    }
    // copy on the components
    for ( UInt icomp = 0; icomp < nb; icomp++ )
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst ( i );
            mat_icomp ( iloc, iloc ) += mat_tmp ( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst ( i );
            jloc = fe.patternSecond ( i );
            mat_icomp ( iloc, jloc ) += mat_tmp ( iloc, jloc );
            mat_icomp ( jloc, iloc ) += mat_tmp ( jloc, iloc );
        }
    }
}

//
// coeff[q]*Mass
//
void mass ( const std::vector<Real>& coef, MatrixElemental& elmat, const CurrentFE& fe,
            int iblock, int jblock, UInt nb )
/*
  Mass matrix: \int v_i v_j (nb blocks on the diagonal, nb>1)
*/
{
    Matrix mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );
    UInt i, ig;
    int iloc, jloc;
    Real s;//, coef_s;
    mat_tmp = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += coef[ig] * fe.phi ( iloc, ig ) * fe.phi ( iloc, ig ) * fe.weightDet ( ig );
        }
        mat_tmp ( iloc, iloc ) = s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += coef[ig] * fe.phi ( iloc, ig ) * fe.phi ( jloc, ig ) * fe.weightDet ( ig );
        }
        mat_tmp ( iloc, jloc ) += s;
        mat_tmp ( jloc, iloc ) += s;
    }
    // copy on the components
    for ( UInt icomp = 0; icomp < nb; icomp++ )
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst ( i );
            mat_icomp ( iloc, iloc ) += mat_tmp ( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst ( i );
            jloc = fe.patternSecond ( i );
            mat_icomp ( iloc, jloc ) += mat_tmp ( iloc, jloc );
            mat_icomp ( jloc, iloc ) += mat_tmp ( jloc, iloc );
        }
    }
}






void stiff_divgrad ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    // div u^k at quad pts
    std::vector<Real> duk (fe.nbQuadPt() );

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        s = 0;
        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {
            // s = 0.0; Alessandro
            for ( UInt i = 0; i < fe.nbFEDof(); i++ )
            {
                s += fe.phiDer ( i, icoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    // costruzione di \div u^k at a quadrature point
            }

            // duk[ ig ] = s; Alessandro
        }// chiude il ciclo su icoor
        duk[ ig ] = s;
    }// chiude il ciclo su ig

    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0.0;
            for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
            {
                for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += duk[ ig ] * fe.phiDer ( i, k, ig ) *  fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }

}

// Stiffness matrix: coef * ( (\div u) \grad u_k : \grad v  ) controllato!!!
void stiff_divgrad_2 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    boost::multi_array<Real, 3> guk (boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {
            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    //  \grad u^k at a quadrature point
                }
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }

    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                        {
                            s += fe.phiDer ( j, jcoor, ig ) * guk[ icoor ][ k ][ ig ] * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                        }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}

// Stiffness matrix: coef * ( \grad u_k : \grad u_k) *( \grad u : \grad v  ) controllato!!!
void stiff_gradgrad ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s, s1;
    std::vector<Real> gguk (fe.nbQuadPt() );

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        s = 0;
        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {
            for ( UInt l = 0; l < fe.nbLocalCoor(); l++ )
            {
                s1 = 0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s1 += fe.phiDer ( i, l, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
                }
                s += s1 * s1;
            }
        }// chiude il ciclo su icoor
        gguk[ ig ] = s;
    }// chiude il ciclo su ig

    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0.0;
            for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
            {
                for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                {
                    s += gguk[ ig ] * fe.phiDer ( i, k, ig ) *  fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }
}

// Stiffness matrix: coef * ( \grad u_k : \grad u) *( \grad u_k : \grad v  ) controllato!!!
void stiff_gradgrad_2 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    boost::multi_array<Real, 3> guk (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {
            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    //  \grad u^k at a quadrature point
                }
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
                    {
                        for ( UInt l = 0; l < fe.nbLocalCoor(); ++l )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += guk[ jcoor ][ l ][ ig ] * fe.phiDer ( j, l, ig ) * guk[ icoor ][ k ][ ig ] * fe.phiDer (i, k, ig ) * fe.weightDet ( ig );    // il ciclo sui nodi di quadratura
                            }
                        }                                                                                                                                                                    // fa l'integrale
                    }
                    mat ( i, j ) += coef  * s;
                }
            }
        }
    }
}

// Stiffness matrix: coef * ( \grad d^k \grad d : \grad v  )controllato!!!
void stiff_dergrad_gradbis ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    boost::multi_array<Real, 3> guk (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    //  \grad u^k at a quadrature point
                }
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
                    {
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                        {
                            s += guk[ icoor ][ jcoor ][ ig ] * fe.phiDer ( i, k, ig ) *  fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}

// Stiffness matrix: coef * ( \grad u^k [\grad d]^T : \grad v  ) controllato!!!
void stiff_dergrad_gradbis_Tr ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    boost::multi_array<Real, 3> guk (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    //  \grad u^k at a quadrature point
                }
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0.0;
                    for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
                    {
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                        {
                            s += guk[ icoor ][ k ][ ig ]  * fe.phiDer ( j, k, ig ) * fe.phiDer ( i, jcoor, ig ) * fe.weightDet ( ig );
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}

// Stiffness matrix: coef * ( \grad \delta d \grad d^k : \grad v  ) controllato!!!
void stiff_dergrad_gradbis_2 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    boost::multi_array<Real, 3> guk (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    //  \grad u^k at a quadrature point
                }
                guk[ icoor ][ jcoor ][ ig ] = s;
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
            for ( UInt l = 0; l < fe.nbLocalCoor(); ++l )
            {
                for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
                {
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    {
                        s += guk[ l ][ k ][ ig ] * fe.phiDer ( i, k, ig ) *  fe.phiDer ( j, l, ig ) * fe.weightDet ( ig );
                    }
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }
}

// Stiffness matrix: coef * ( \grad \delta d [\grad d^k]^T : \grad v  )
void stiff_dergrad_gradbis_Tr_2 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    boost::multi_array<Real, 3> guk (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    //  \grad u^k at a quadrature point
                }
                guk[ icoor ][ jcoor ][ ig ] = s;
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
            for ( UInt l = 0; l < fe.nbLocalCoor(); ++l )
            {
                for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
                {
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    {
                        s += guk[ k ][ l ][ ig ] * fe.phiDer ( i, k, ig ) *  fe.phiDer ( j, l, ig ) * fe.weightDet ( ig );
                    }
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor ) // copia del blocco sulla diagonale
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }
}

//  Stiffness matrix: coef * (  \grad u^k [\grad u^k]^T \grad u : \grad v  ) controllato!!!
void stiff_gradgradTr_gradbis ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    boost::multi_array<Real, 3> guk_gukT (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt n = 0; n < fe.nbLocalCoor(); n++ )
                {
                    for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    {
                        for ( UInt j = 0; j < fe.nbFEDof(); j++ )
                        {
                            s  += fe.phiDer ( i, n, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ] * fe.phiDer ( j, n, ig ) * uk_loc.vec() [ j + jcoor * fe.nbFEDof() ] ;    // \grad u^k  [\grad u^k]^T  at each quadrature point
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

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor ); // estrae il blocco (icoor, jcoor)

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                        {
                            s += fe.phiDer ( i, k, ig ) * guk_gukT[ icoor ][ jcoor ][ ig ] * fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                        }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}

//  Stiffness matrix: coef * (  \grad u^k [\grad u]^T \grad u^k : \grad v  )controllato!!!
void stiff_gradgradTr_gradbis_2 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    boost::multi_array<Real, 3> guk (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    //  \grad u^k at a quadrature point
                }
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }

    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor ); // estrae il blocco (icoor, jcoor)

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt l = 0; l < fe.nbLocalCoor(); ++l )
                    {
                        for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
                        {
                            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                            {
                                s += guk[ icoor ][ l ][ ig ] * guk[ jcoor ][ k ][ ig ] * fe.phiDer ( i, k, ig ) *  fe.phiDer ( j, l, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}

//  Stiffness matrix: coef * (  \grad u [\grad u^k]^T \grad u^k : \grad v  )controllato!!!
void stiff_gradgradTr_gradbis_3 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    boost::multi_array<Real, 3> guk_gukT (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    // attenzione in questa funzione si deve usare il trasposto
    // loop on quadrature points                                                // (\grad u^k  [\grad u^k]^T )^T
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt n = 0; n < fe.nbLocalCoor(); n++ )
                {
                    for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                    {
                        for ( UInt j = 0; j < fe.nbFEDof(); j++ )
                        {
                            s  += fe.phiDer ( i, n, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ] * fe.phiDer ( j, n, ig ) * uk_loc.vec() [ j + jcoor * fe.nbFEDof() ] ;    // \grad u^k  [\grad u^k]^T  at each quadrature point
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
            for ( UInt l = 0; l < fe.nbLocalCoor(); ++l )
            {
                for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
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

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor ) // copia del blocco sulla diagonale
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }
}







void ipstab_grad ( const Real         coef,
                   MatrixElemental&           elmat,
                   const CurrentFE&   fe1,
                   const CurrentFE&   fe2,
                   const CurrentFEManifold& bdfe,
                   int iblock, int jblock )
{
    /*
      Interior penalty stabilization: coef*\int_{face} grad u1_i . grad v1_j
    */

    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );



    Real sum, sum1, sum2;
    UInt icoor, jcoor, i, j;
    UInt ig;
    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];

    boost::multi_array<Real, 3> phid1 (
        boost::extents[fe1.nbFEDof()][fe1.nbLocalCoor()][bdfe.nbQuadPt()]);
    boost::multi_array<Real, 3> phid2 (
        boost::extents[fe2.nbFEDof()][fe2.nbLocalCoor()][bdfe.nbQuadPt()]);

    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap ( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap ( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {
        // first derivatives on quadrature points
        bdfe.coorQuadPt ( x[ 0 ], x[ 1 ], x[ 2 ], ig );      // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi ( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi ( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac ( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac ( icoor, jcoor, 0 ) * drp2[ jcoor ];
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
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt() ; ++ig )
                {
                    sum += phid1[ i ][ icoor ][ ig ] * phid2[ j ][ icoor ][ ig ] * bdfe.wRootDetMetric ( ig );
                }
            mat ( i, j ) += coef * sum;
        }
    }

}






void ipstab_grad ( const Real         coef,
                   MatrixElemental&           elmat,
                   const CurrentFE&   fe1,
                   const CurrentFE&   fe2,
                   const CurrentFEManifold& bdfe,
                   int iblock, int jblock,
                   int nb )
{
    /*
      Interior penalty stabilization: coef*\int_{face} grad u1_i . grad v1_j
    */


    MatrixElemental::matrix_type mat_tmp ( fe1.nbFEDof(), fe2.nbFEDof() );


    Real sum, sum1, sum2;
    UInt icoor, jcoor, i, j;
    UInt ig;
    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];

    boost::multi_array<Real, 3> phid1 (
        boost::extents[fe1.nbFEDof()][fe1.nbLocalCoor()][bdfe.nbQuadPt()]);
    boost::multi_array<Real, 3> phid2 (
        boost::extents[fe2.nbFEDof()][fe2.nbLocalCoor()][bdfe.nbQuadPt()]);

    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap ( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap ( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {
        // first derivatives on quadrature points
        bdfe.coorQuadPt ( x[ 0 ], x[ 1 ], x[ 2 ], ig );      // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi ( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi ( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac ( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac ( icoor, jcoor, 0 ) * drp2[ jcoor ];
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
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt() ; ++ig )
                {
                    sum += phid1[ i ][ icoor ][ ig ] * phid2[ j ][ icoor ][ ig ] * bdfe.wRootDetMetric ( ig );
                }
            mat_tmp ( i, j ) = coef * sum;
        }
    }


    // copy on the components
    for ( int icomp = 0; icomp < nb; icomp++ )
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icomp, jblock + icomp );
        mat_icomp += mat_tmp;
    }
}





void ipstab_bgrad ( const Real         coef,
                    MatrixElemental&           elmat,
                    const CurrentFE&   fe1,
                    const CurrentFE&   fe2,
                    const VectorElemental&     beta,
                    const CurrentFEManifold& bdfe,
                    int iblock, int jblock,
                    int nb )
{
    /*
      Interior penalty stabilization: coef*\int_{face} (\beta1 . grad u1_i) . (\beta2 . grad v2_j)
    */

    MatrixElemental::matrix_type mat_tmp ( fe1.nbFEDof(), fe2.nbFEDof() );

    Real sum, sum1, sum2;
    UInt i, j;
    UInt icoor, jcoor;
    UInt ig;

    //
    // convection velocity \beta on the boundary quadrature points
    //
    boost::multi_array<Real, 2> b (
        boost::extents[fe1.nbLocalCoor()][bdfe.nbQuadPt()]);

    for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
    {
        for ( ig = 0; ig < bdfe.nbQuadPt(); ig++ )
        {
            sum = 0;
            for ( i = 0; i < bdfe.nbFEDof(); ++i )
            {
                sum += bdfe.phi ( i, ig ) * beta.vec() [ icoor * bdfe.nbLocalCoor() + i ];
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

    boost::multi_array<Real, 3> phid1 (
        boost::extents[fe1.nbFEDof()][fe1.nbLocalCoor()][bdfe.nbQuadPt()]);
    boost::multi_array<Real, 3> phid2 (
        boost::extents[fe2.nbFEDof()][fe2.nbLocalCoor()][bdfe.nbQuadPt()]);

    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap ( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap ( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {
        // first derivatives on quadrature points
        bdfe.coorQuadPt ( x[ 0 ], x[ 1 ], x[ 2 ], ig );      // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi ( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi ( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac ( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac ( icoor, jcoor, 0 ) * drp2[ jcoor ];
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
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
                for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
                    for ( ig = 0; ig < bdfe.nbQuadPt(); ig++ )
                        sum += phid1[ i ][ icoor ][ ig ] * phid2[ j ][ jcoor ][ ig ]
                               * b[ icoor ][ ig ] * b[ jcoor ][ ig ]
                               * bdfe.wRootDetMetric ( ig );
            mat_tmp ( i, j ) = coef * sum;
        }
    }

    // copy on the components
    for ( int icomp = 0; icomp < nb; icomp++ )
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icomp, jblock + icomp );
        mat_icomp += mat_tmp;
    }

}





void ipstab_div ( const Real coef, MatrixElemental& elmat, const CurrentFE& fe1, const CurrentFE& fe2,
                  const CurrentFEManifold& bdfe, int iblock, int jblock )
{
    /*
      Interior penalty stabilization: coef*\int_{face} div u . div v
    */

    Real sum, sum1, sum2;
    UInt i, j, icoor, jcoor;
    UInt ig;
    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];

    boost::multi_array<Real, 3> phid1 (
        boost::extents[fe1.nbFEDof()][fe1.nbLocalCoor()][bdfe.nbQuadPt()]);
    boost::multi_array<Real, 3> phid2 (
        boost::extents[fe2.nbFEDof()][fe2.nbLocalCoor()][bdfe.nbQuadPt()]);

    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap ( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap ( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {
        // first derivatives on quadrature points
        bdfe.coorQuadPt ( x[ 0 ], x[ 1 ], x[ 2 ], ig );      // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi ( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi ( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac ( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac ( icoor, jcoor, 0 ) * drp2[ jcoor ];
                }
                phid1[ i ][ icoor ][ ig ] = sum1;
                phid2[ i ][ icoor ][ ig ] = sum2;
            }
        }
    }

    for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
    {
        for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
        {
            MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icoor, jblock + jcoor );
            // Loop on rows
            for ( i = 0; i < fe1.nbFEDof(); ++i )
            {
                // Loop on columns
                for ( j = 0; j < fe2.nbFEDof(); ++j )
                {
                    sum = 0.0;
                    for ( ig = 0; ig < bdfe.nbQuadPt(); ig++ )
                    {
                        sum += phid1[ i ][ icoor ][ ig ] * phid2[ j ][ jcoor ][ ig ] * bdfe.wRootDetMetric ( ig );
                    }
                    mat_icomp ( i, j ) += coef * sum;
                }
            }
        }
    }

}

void ipstab_bagrad ( const Real coef, MatrixElemental& elmat,
                     const CurrentFE& fe1, const CurrentFE& fe2,
                     const VectorElemental& beta, const CurrentFEManifold& bdfe,
                     int iblock, int jblock )
{

    // Interior penalty stabilization:
    // coef < |\beta . n|^2 / |\beta| \grad p1, \grad q2 >

    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );

    Real sum, sum1, sum2;
    UInt icoor, jcoor;
    UInt i, j, ig;

    boost::multi_array<Real, 3> phid1 (
        boost::extents[fe1.nbFEDof()][fe1.nbLocalCoor()][bdfe.nbQuadPt()]);
    boost::multi_array<Real, 3> phid2 (
        boost::extents[fe2.nbFEDof()][fe2.nbLocalCoor()][bdfe.nbQuadPt()]);

    std::vector<Real> x (3), rx1 (3), drp1 (3), rx2 (3), drp2 (3);
    std::vector<Real> b1 (3), b2 (3);

    fe1.coorMap ( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap ( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    //
    // convection velocity term |\beta . n|^2 / |\beta|
    // on the boundary quadrature points
    //
    std::vector<Real> ba2 (bdfe.nbQuadPt() );

    for ( ig = 0; ig < bdfe.nbQuadPt(); ig++ )
    {
        sum1 = 0;
        sum2 = 0;
        for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
        {
            for ( i = 0; i < bdfe.nbLocalCoor(); ++i )
            {
                Real betaLoc = bdfe.phi ( i, ig ) *
                               beta.vec() [ icoor * bdfe.nbLocalCoor() + i ];
                sum1 += betaLoc * bdfe.normal (icoor, ig);
                sum2 += betaLoc * betaLoc;
            }
        }
        ba2[ ig ] = sum2 == 0 ? 0 : sum1 * sum1 / std::pow ( sum2, 0.5 );
    }

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {
        // first derivatives on quadrature points
        bdfe.coorQuadPt ( x[ 0 ], x[ 1 ], x[ 2 ], ig );      // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi ( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi ( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac ( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac ( icoor, jcoor, 0 ) * drp2[ jcoor ];
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
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt() ; ++ig )
                    sum += ba2[ ig ] *
                           phid1[ i ][ icoor ][ ig ] *
                           phid2[ j ][ icoor ][ ig ] *
                           bdfe.wRootDetMetric ( ig );
            mat ( i, j ) += coef * sum;
        }
    }

}


// coef < |\beta . n| \grad p1, \grad q2 >
// p1   lives in fe1
// q2   lives in fe2
// beta lives in fe3


void ipstab_bagrad ( const Real         coef,
                     MatrixElemental&           elmat,
                     const CurrentFE&   fe1,
                     const CurrentFE&   fe2,
                     const CurrentFE&   fe3,
                     const VectorElemental&     beta,
                     const CurrentFEManifold& bdfe,
                     int iblock, int    jblock )
{

    // Interior penalty stabilization:
    // coef < |\beta.n| \grad p1, \grad q2 >

    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );

    Real sum, sum1, sum2;
    UInt icoor, jcoor;
    UInt i, j, ig;

    boost::multi_array<Real, 3> phid1 (
        boost::extents[fe1.nbFEDof()][fe1.nbLocalCoor()][bdfe.nbQuadPt()]);
    boost::multi_array<Real, 3> phid2 (
        boost::extents[fe2.nbFEDof()][fe2.nbLocalCoor()][bdfe.nbQuadPt()]);

    std::vector<Real> x (3), rx1 (3), drp1 (3), rx2 (3), drp2 (3);
    std::vector<Real> b1 (3), b2 (3);

    fe1.coorMap ( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap ( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    //
    // convection velocity term |\beta . n|
    // on the boundary quadrature points
    //
    std::vector<Real> bn (bdfe.nbQuadPt() );

    for ( ig = 0; ig < bdfe.nbQuadPt(); ig++ )
    {
        sum1 = 0;
        sum2 = 0;

        for ( icoor = 0; icoor < fe3.nbLocalCoor(); ++icoor )
        {
            for ( i = 0; i < fe3.nbFEDof(); ++i )
            {
                Real betaLoc = fe3.phi ( i, ig ) *
                               beta.vec() [ icoor * fe3.nbFEDof() + i ];
                sum1 += betaLoc * bdfe.normal (icoor, ig);
            }
        }
        bn[ ig ] = std::fabs (sum1);
    }

    for ( UInt ig = 0; ig < bdfe.nbQuadPt(); ++ig )
    {
        // first derivatives on quadrature points
        bdfe.coorQuadPt ( x[ 0 ], x[ 1 ], x[ 2 ], ig );      // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
            {
                sum1 += fe1.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac ( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbFEDof(); ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                drp1[ icoor ] = fe1.refFE().dPhi ( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE().dPhi ( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbLocalCoor(); ++jcoor )
                {
                    sum1 += fe1.tInvJac ( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac ( icoor, jcoor, 0 ) * drp2[ jcoor ];
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
            for ( icoor = 0; icoor < fe1.nbLocalCoor(); ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt() ; ++ig )
                    sum += bn[ ig ] *
                           phid1[ i ][ icoor ][ ig ] *
                           phid2[ j ][ icoor ][ ig ] *
                           bdfe.wRootDetMetric ( ig );
            mat ( i, j ) += coef * sum;
        }
    }

}


void stiff ( Real coef, MatrixElemental& elmat, const CurrentFE& fe,
             int iblock, int jblock )
/*
  Stiffness matrix: coef*\int grad v_i . grad v_j
*/
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    UInt iloc, jloc;
    UInt i, icoor, ig;
    double s, coef_s;
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                s += fe.phiDer ( iloc, icoor, ig ) * fe.phiDer ( iloc, icoor, ig )
                     * fe.weightDet ( ig );
        }
        mat ( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                s += fe.phiDer ( iloc, icoor, ig ) * fe.phiDer ( jloc, icoor, ig ) *
                     fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat ( iloc, jloc ) += coef_s;
        mat ( jloc, iloc ) += coef_s;
    }
}


void stiff ( Real coef, Real ( *fct ) ( Real, Real, Real ), MatrixElemental& elmat,
             const CurrentFE& fe, int iblock, int jblock )
/*
  Stiffness matrix: coef*\int grad v_i . grad v_j
*/
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    UInt iloc, jloc;
    UInt i, icoor, ig;
    double s, coef_s, coef_f;
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            coef_f = fct ( fe.quadPt ( ig, 0 ), fe.quadPt ( ig, 1 ), fe.quadPt ( ig, 2 ) );
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                s += coef_f * fe.phiDer ( iloc, icoor, ig ) * fe.phiDer ( iloc, icoor, ig )
                     * fe.weightDet ( ig );
        }
        mat ( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            coef_f = fct ( fe.quadPt ( ig, 0 ), fe.quadPt ( ig, 1 ), fe.quadPt ( ig, 2 ) );
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                s += coef_f * fe.phiDer ( iloc, icoor, ig ) * fe.phiDer ( jloc, icoor, ig ) *
                     fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat ( iloc, jloc ) += coef_s;
        mat ( jloc, iloc ) += coef_s;
    }
}
//

void stiff ( Real coef, MatrixElemental& elmat, const CurrentFE& fe,
             int iblock, int jblock, int nb )
/*
  Stiffness matrix: coef*\int grad v_i . grad v_j (nb blocks on the diagonal, nb>1)
*/
{


    Matrix mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );

    UInt iloc, jloc;
    UInt i, icoor, ig;
    double s, coef_s;
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                s += fe.phiDer ( iloc, icoor, ig ) * fe.phiDer ( iloc, icoor, ig )
                     * fe.weightDet ( ig );
        }
        mat_tmp ( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                s += fe.phiDer ( iloc, icoor, ig ) * fe.phiDer ( jloc, icoor, ig ) *
                     fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat_tmp ( iloc, jloc ) += coef_s;
        mat_tmp ( jloc, iloc ) += coef_s;
    }
    // copy on the components
    for ( int icomp = 0; icomp < nb; icomp++ )
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst ( i );
            mat_icomp ( iloc, iloc ) += mat_tmp ( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst ( i );
            jloc = fe.patternSecond ( i );
            mat_icomp ( iloc, jloc ) += mat_tmp ( iloc, jloc );
            mat_icomp ( jloc, iloc ) += mat_tmp ( jloc, iloc );
        }
    }
}


void stiff ( const std::vector<Real>& coef, MatrixElemental& elmat, const CurrentFE& fe,
             int iblock, int jblock, int nb )
/*
  Stiffness matrix: coef*\int grad v_i . grad v_j (nb blocks on the diagonal, nb>1)
  with coef given in a vector (one element per quadrature point)
*/
{

    Matrix mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );

    UInt iloc, jloc;
    UInt i, icoor, ig;
    double s;//, coef_s;
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                s += coef[ig] * fe.phiDer ( iloc, icoor, ig ) * fe.phiDer ( iloc, icoor, ig )
                     * fe.weightDet ( ig );
        }
        mat_tmp ( iloc, iloc ) += s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                s += coef[ig] * fe.phiDer ( iloc, icoor, ig ) * fe.phiDer ( jloc, icoor, ig ) *
                     fe.weightDet ( ig );
        }
        mat_tmp ( iloc, jloc ) += s;
        mat_tmp ( jloc, iloc ) += s;
    }
    // copy on the components
    for ( int icomp = 0; icomp < nb; icomp++ )
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst ( i );
            mat_icomp ( iloc, iloc ) += mat_tmp ( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst ( i );
            jloc = fe.patternSecond ( i );
            mat_icomp ( iloc, jloc ) += mat_tmp ( iloc, jloc );
            mat_icomp ( jloc, iloc ) += mat_tmp ( jloc, iloc );
        }
    }
}


// VC - December 2004
//
void stiff_curl ( Real coef, MatrixElemental& elmat, const CurrentFE& fe,
                  int iblock, int jblock, int /*nb*/ )


/*
  Stiffness matrix: coef*\int curl v_i . curl v_j (nb blocks on the diagonal, nb>1)
*/
{


    Matrix mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp11 ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp11 = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp12 ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp12 = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp13 ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp13 = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp21 ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp21 = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp22 ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp22 = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp23 ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp23 = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp31 ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp31 = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp32 ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp32 = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );
    Matrix mat_tmp33 ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp33 = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );



    UInt iloc, jloc;
    UInt i, ig;
    double s, coef_s;

    // diagonal 11
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer ( iloc, 1, ig ) * fe.phiDer ( iloc, 1, ig ) * fe.weightDet ( ig )
                + fe.phiDer ( iloc, 2, ig ) * fe.phiDer ( iloc, 2, ig ) * fe.weightDet ( ig ) ;
        }
        mat_tmp11 ( iloc, iloc ) += coef * s;
    }
    // extra diagonal 11
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer ( iloc, 1, ig ) * fe.phiDer ( jloc, 1, ig ) * fe.weightDet ( ig )
                + fe.phiDer ( iloc, 2, ig ) * fe.phiDer ( jloc, 2, ig ) * fe.weightDet ( ig )       ;
        }
        coef_s = coef * s;
        mat_tmp11 ( iloc, jloc ) += coef_s;
        mat_tmp11 ( jloc, iloc ) += coef_s;
    }

    // diagonal 12
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 1, ig ) * fe.phiDer ( iloc, 0, ig ) * fe.weightDet ( ig );
        }
        mat_tmp12 ( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 12
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 1, ig ) * fe.phiDer ( jloc, 0, ig ) * fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat_tmp12 ( iloc, jloc ) -= coef_s;
        mat_tmp12 ( jloc, iloc ) -= coef_s;
    }

    // diagonal 13
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 2, ig ) * fe.phiDer ( iloc, 0, ig ) * fe.weightDet ( ig );
        }
        mat_tmp13 ( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 13
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 2, ig ) * fe.phiDer ( jloc, 0, ig ) * fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat_tmp13 ( iloc, jloc ) -= coef_s;
        mat_tmp13 ( jloc, iloc ) -= coef_s;
    }

    // diagonal 21
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 0, ig ) * fe.phiDer ( iloc, 1, ig ) * fe.weightDet ( ig );
        }
        mat_tmp21 ( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 21
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 0, ig ) * fe.phiDer ( jloc, 1, ig ) * fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat_tmp21 ( iloc, jloc ) -= coef_s;
        mat_tmp21 ( jloc, iloc ) -= coef_s;
    }

    // diagonal 22
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer ( iloc, 0, ig ) * fe.phiDer ( iloc, 0, ig ) * fe.weightDet ( ig )
                + fe.phiDer ( iloc, 2, ig ) * fe.phiDer ( iloc, 2, ig ) * fe.weightDet ( ig ) ;
        }
        mat_tmp22 ( iloc, iloc ) += coef * s;
    }
    // extra diagonal 22
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer ( iloc, 0, ig ) * fe.phiDer ( jloc, 0, ig ) * fe.weightDet ( ig )
                + fe.phiDer ( iloc, 2, ig ) * fe.phiDer ( jloc, 2, ig ) * fe.weightDet ( ig )       ;
        }
        coef_s = coef * s;
        mat_tmp22 ( iloc, jloc ) += coef_s;
        mat_tmp22 ( jloc, iloc ) += coef_s;
    }

    // diagonal 23
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 2, ig ) * fe.phiDer ( iloc, 1, ig ) * fe.weightDet ( ig );
        }
        mat_tmp23 ( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 23
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 2, ig ) * fe.phiDer ( jloc, 1, ig ) * fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat_tmp23 ( iloc, jloc ) -= coef_s;
        mat_tmp23 ( jloc, iloc ) -= coef_s;
    }

    // diagonal 31
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 0, ig ) * fe.phiDer ( iloc, 2, ig ) * fe.weightDet ( ig );
        }
        mat_tmp31 ( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 31
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 0, ig ) * fe.phiDer ( jloc, 2, ig ) * fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat_tmp31 ( iloc, jloc ) -= coef_s;
        mat_tmp31 ( jloc, iloc ) -= coef_s;
    }

    // diagonal 32
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 1, ig ) * fe.phiDer ( iloc, 2, ig ) * fe.weightDet ( ig );
        }
        mat_tmp32 ( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 32
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = - fe.phiDer ( iloc, 1, ig ) * fe.phiDer ( jloc, 2, ig ) * fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat_tmp32 ( iloc, jloc ) -= coef_s;
        mat_tmp32 ( jloc, iloc ) -= coef_s;
    }

    // diagonal 33
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer ( iloc, 0, ig ) * fe.phiDer ( iloc, 0, ig ) * fe.weightDet ( ig )
                + fe.phiDer ( iloc, 1, ig ) * fe.phiDer ( iloc, 1, ig ) * fe.weightDet ( ig ) ;
        }
        mat_tmp33 ( iloc, iloc ) += coef * s;
    }
    // extra diagonal 33
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s = fe.phiDer ( iloc, 1, ig ) * fe.phiDer ( jloc, 1, ig ) * fe.weightDet ( ig )
                + fe.phiDer ( iloc, 2, ig ) * fe.phiDer ( jloc, 2, ig ) * fe.weightDet ( ig )       ;
        }
        coef_s = coef * s;
        mat_tmp33 ( iloc, jloc ) += coef_s;
        mat_tmp33 ( jloc, iloc ) += coef_s;
    }

    MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + 0, jblock + 0 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        mat_icomp ( iloc, iloc ) += mat_tmp11 ( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        mat_icomp ( iloc, jloc ) += mat_tmp11 ( iloc, jloc );
        mat_icomp ( jloc, iloc ) += mat_tmp11 ( jloc, iloc );
    }

    mat_icomp = elmat.block ( iblock + 0, jblock + 1 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        mat_icomp ( iloc, iloc ) -= mat_tmp12 ( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        mat_icomp ( iloc, jloc ) -= mat_tmp12 ( iloc, jloc );
        mat_icomp ( jloc, iloc ) -= mat_tmp12 ( jloc, iloc );
    }

    mat_icomp = elmat.block ( iblock + 0, jblock + 2 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        mat_icomp ( iloc, iloc ) -= mat_tmp13 ( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        mat_icomp ( iloc, jloc ) -= mat_tmp13 ( iloc, jloc );
        mat_icomp ( jloc, iloc ) -= mat_tmp13 ( jloc, iloc );
    }

    mat_icomp = elmat.block ( iblock + 1, jblock + 0 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        mat_icomp ( iloc, iloc ) -= mat_tmp21 ( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        mat_icomp ( iloc, jloc ) -= mat_tmp21 ( iloc, jloc );
        mat_icomp ( jloc, iloc ) -= mat_tmp21 ( jloc, iloc );
    }

    mat_icomp = elmat.block ( iblock + 1, jblock + 1 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        mat_icomp ( iloc, iloc ) += mat_tmp22 ( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        mat_icomp ( iloc, jloc ) += mat_tmp22 ( iloc, jloc );
        mat_icomp ( jloc, iloc ) += mat_tmp22 ( jloc, iloc );
    }

    mat_icomp = elmat.block ( iblock + 1, jblock + 2 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        mat_icomp ( iloc, iloc ) -= mat_tmp23 ( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        mat_icomp ( iloc, jloc ) -= mat_tmp23 ( iloc, jloc );
        mat_icomp ( jloc, iloc ) -= mat_tmp23 ( jloc, iloc );
    }

    mat_icomp = elmat.block ( iblock + 2, jblock + 0 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        mat_icomp ( iloc, iloc ) -= mat_tmp31 ( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        mat_icomp ( iloc, jloc ) -= mat_tmp31 ( iloc, jloc );
        mat_icomp ( jloc, iloc ) -= mat_tmp31 ( jloc, iloc );
    }

    mat_icomp = elmat.block ( iblock + 2, jblock + 1 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        mat_icomp ( iloc, iloc ) -= mat_tmp32 ( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        mat_icomp ( iloc, jloc ) -= mat_tmp32 ( iloc, jloc );
        mat_icomp ( jloc, iloc ) -= mat_tmp32 ( jloc, iloc );
    }

    mat_icomp = elmat.block ( iblock + 2, jblock + 2 );
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        mat_icomp ( iloc, iloc ) += mat_tmp33 ( iloc, iloc );
    }
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        mat_icomp ( iloc, jloc ) += mat_tmp33 ( iloc, jloc );
        mat_icomp ( jloc, iloc ) += mat_tmp ( jloc, iloc );
    }
}


/*
  Stiffness matrix: coef * ( div u , div v )
*/
void stiff_div ( Real coef, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;

    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    {
                        s += fe.phiDer ( i, icoor, ig ) * fe.phiDer ( j, jcoor, ig ) * fe.weightDet ( ig );
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}



/*
  Stiffness matrix: coef * ( [\grad u^k]^T \grad d : \grad v  )
*/
void stiff_dergradbis ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;
    boost::multi_array<Real, 3> guk (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    //  \grad u^k at a quadrature point
                }
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                        {
                            s += fe.phiDer ( i, k, ig ) * guk[ jcoor ][ icoor ][ ig ] * fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                        }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}




/*
  Stiffness matrix: coef * ( [\grad u]^T \grad u^k [\grad u^k]^T \grad u : \grad v  ) for Newton on St-Venant
*/
void stiff_dergrad ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{

    double s;

    boost::multi_array<Real, 3> guk (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    //  \grad u^k at a quadrature point
                }
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                        {
                            s += fe.phiDer ( i, k, ig ) * ( guk[ jcoor ][ k ][ ig ] * fe.phiDer ( j, icoor, ig )
                                                            + guk[ jcoor ][ icoor ][ ig ] * fe.phiDer ( j, k, ig ) ) * fe.weightDet ( ig );
                        }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}





//
// coef * ( \tr { [\grad u^k]^T \grad u }, \div v  ) for Newton on St-Venant
//
//
void stiff_derdiv ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe )
{
    boost::multi_array<Real, 3> guk (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbQuadPt()]);

    Real s;

    // loop on quadrature points
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < fe.nbFEDof(); i++ )
                {
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    //  \grad u^k at a quadrature point
                }
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt k = 0; k < fe.nbLocalCoor(); ++k )
                        for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
                        {
                            s += fe.phiDer ( i, icoor, ig ) * guk[ jcoor ][ k ][ ig ] * fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                        }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}




void stiff_strain ( Real coef, MatrixElemental& elmat, const CurrentFE& fe )
/*
  Stiffness matrix: coef * ( e(u) , e(v) )
*/
{
    double s;
    double tmp = coef * 0.5;

    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0;
            for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
                {
                    s += fe.phiDer ( i, icoor, ig ) * fe.phiDer ( j, icoor, ig ) * fe.weightDet ( ig );
                }
            mat_tmp ( i, j ) = tmp * s;
        }
    }
    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }

    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( UInt jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( UInt j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    {
                        s += fe.phiDer ( i, jcoor, ig ) * fe.phiDer ( j, icoor, ig ) * fe.weightDet ( ig );
                    }
                    mat ( i, j ) += tmp * s;
                }
            }
        }
    }
}

void mass_divw ( Real coef, const VectorElemental& w_loc, MatrixElemental& elmat, const CurrentFE& fe,
                 int iblock, int jblock, UInt nb )
/*
  modified mass matrix: ( div w u,v )
*/
{

    Matrix mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );

    UInt i, icomp, ig, icoor, iloc, jloc;
    Real s, coef_s;
    std::vector<Real> divw (fe.nbQuadPt() );

    // divw at quadrature nodes
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        divw[ ig ] = 0.0;
        for ( icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
            for ( i = 0; i < fe.nbFEDof(); ++i )
            {
                divw[ ig ] += fe.phiDer ( i, icoor, ig ) * w_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
    }

    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += divw[ ig ] * fe.phi ( iloc, ig ) * fe.phi ( iloc, ig ) * fe.weightDet ( ig );
        }
        mat_tmp ( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += divw[ ig ] * fe.phi ( iloc, ig ) * fe.phi ( jloc, ig ) * fe.weightDet ( ig );
        }
        coef_s = coef * s;
        mat_tmp ( iloc, jloc ) += coef_s;
        mat_tmp ( jloc, iloc ) += coef_s;
    }
    // copy on the components
    for ( icomp = 0; icomp < nb; icomp++ )
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst ( i );
            mat_icomp ( iloc, iloc ) += mat_tmp ( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst ( i );
            jloc = fe.patternSecond ( i );
            mat_icomp ( iloc, jloc ) += mat_tmp ( iloc, jloc );
            mat_icomp ( jloc, iloc ) += mat_tmp ( jloc, iloc );
        }
    }
}


void mass_divw (const std::vector<Real>& coef, const VectorElemental& w_loc, MatrixElemental& elmat, const CurrentFE& fe,
                int iblock, int jblock, UInt nb )
/*
  modified mass matrix: ( div w u,v )
*/
{

    Matrix mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp = ZeroMatrix ( fe.nbFEDof(), fe.nbFEDof() );

    UInt i, icomp, ig, icoor, iloc, jloc;
    Real s;
    std::vector<Real> divw (fe.nbQuadPt() );

    // divw at quadrature nodes
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        divw[ ig ] = 0.0;
        for ( icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
            for ( i = 0; i < fe.nbFEDof(); ++i )
            {
                divw[ ig ] += fe.phiDer ( i, icoor, ig ) * w_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
    }

    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += coef[ig] * divw[ ig ] * fe.phi ( iloc, ig ) * fe.phi ( iloc, ig ) * fe.weightDet ( ig );
        }
        mat_tmp ( iloc, iloc ) += s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += coef[ig] * divw[ ig ] * fe.phi ( iloc, ig ) * fe.phi ( jloc, ig ) * fe.weightDet ( ig );
        }
        mat_tmp ( iloc, jloc ) += s;
        mat_tmp ( jloc, iloc ) += s;
    }
    // copy on the components
    for ( icomp = 0; icomp < nb; icomp++ )
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst ( i );
            mat_icomp ( iloc, iloc ) += mat_tmp ( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst ( i );
            jloc = fe.patternSecond ( i );
            mat_icomp ( iloc, jloc ) += mat_tmp ( iloc, jloc );
            mat_icomp ( jloc, iloc ) += mat_tmp ( jloc, iloc );
        }
    }
}




void mass_gradu ( Real coef, const VectorElemental& u0_loc, MatrixElemental& elmat, const CurrentFE& fe )
/*
  modified mass matrix: ( grad u0 u,v )
*/
{

    UInt ig, icoor, jcoor, i, j;
    Real s;

    boost::multi_array<Real, 3> gu0 (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()][fe.nbLocalCoor()]);

    //
    // grad u0 at quadrature nodes
    //
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        for ( icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
            {
                gu0[ ig ][ icoor ][ jcoor ] = 0.0;
                for ( i = 0; i < fe.nbFEDof(); ++i )
                {
                    gu0[ ig ][ icoor ][ jcoor ] += fe.phiDer ( i, jcoor, ig ) * u0_loc.vec() [ i + icoor * fe.nbFEDof() ];
                }
            }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( icoor = 0; icoor < fe.nbLocalCoor(); ++icoor )
    {
        for ( jcoor = 0; jcoor < fe.nbLocalCoor(); ++jcoor )
        {

            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );

            for ( i = 0; i < fe.nbFEDof(); ++i )
            {
                for ( j = 0; j < fe.nbFEDof(); ++j )
                {
                    s = 0;
                    for ( ig = 0; ig < fe.nbQuadPt(); ++ig )
                    {
                        s += gu0[ ig ][ icoor ][ jcoor ] * fe.phi ( i, ig ) * fe.phi ( j, ig ) * fe.weightDet ( ig );
                    }
                    mat ( i, j ) += coef * s;
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
void stiff_sd ( Real coef, const VectorElemental& vec_loc, MatrixElemental& elmat, const CurrentFE& fe, const CurrentFE& fe2,
                int iblock, int jblock, int nb )
/*
  Stiffness matrix for SD Stabilization
*/
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    UInt iloc, jloc;
    UInt i, icoor, ig, jcoor;
    std::vector<Real> coef_v (fe.nbLocalCoor() );

    double s, coef_s;
    //    int nbN1=fe.nbFEDof();
    UInt nbN2 = fe2.nbFEDof();
    //
    // diagonal
    //
    for ( i = 0; i < fe.nbDiag(); i++ )
    {
        iloc = fe.patternFirst ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            {
                coef_v[ icoor ] = 0.;
            }

            // computation of the convection term in the quadrature nodes
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            {
                for ( UInt iloc = 0; iloc < nbN2; iloc++ )
                {
                    coef_v[ icoor ] += vec_loc.vec() [ iloc + icoor * nbN2 ] * fe2.phi ( iloc, ig );
                }
            }

            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            {
                for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                {
                    s += coef_v[ icoor ] * fe.phiDer ( iloc, icoor, ig ) * coef_v[ jcoor ] * fe.phiDer ( iloc, jcoor, ig )
                         * fe.weightDet ( ig );
                }
            }
        }
        mat ( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
    {
        iloc = fe.patternFirst ( i );
        jloc = fe.patternSecond ( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            {
                coef_v[ icoor ] = 0.;
            }

            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            {
                for ( UInt iloc = 0; iloc < nbN2; iloc++ )
                {
                    coef_v[ icoor ] += vec_loc.vec() [ iloc + icoor * nbN2 ] * fe2.phi ( iloc, ig );
                }
            }

            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            {
                for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                {
                    s += coef_v[ icoor ] * fe.phiDer ( iloc, icoor, ig ) * coef_v[ jcoor ] * fe.phiDer ( jloc, jcoor, ig )
                         * fe.weightDet ( ig );
                }
            }
        }
        coef_s = coef * s;
        mat ( iloc, jloc ) += coef_s;
        mat ( jloc, iloc ) += coef_s;
    }
    // copy on the other components (if necessary, i.e. if nb>1)
    for ( int icomp = 1; icomp < nb; icomp++ )
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icomp, jblock + icomp );
        for ( i = 0; i < fe.nbDiag(); i++ )
        {
            iloc = fe.patternFirst ( i );
            mat_icomp ( iloc, iloc ) += mat ( iloc, iloc );
        }
        for ( i = fe.nbDiag(); i < fe.nbDiag() + fe.nbUpper(); i++ )
        {
            iloc = fe.patternFirst ( i );
            jloc = fe.patternSecond ( i );
            mat_icomp ( iloc, jloc ) += mat ( iloc, jloc );
            mat_icomp ( jloc, iloc ) += mat ( jloc, iloc );
        }
    }
}


//
void grad ( const int icoor, Real coef, MatrixElemental& elmat,
            const CurrentFE& fe_u, const CurrentFE& fe_p,
            int iblock, int jblock )
/*
  \int q_j \frac{\partial v_i}{\partial x_icoor}
*/
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
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

                s -= fe_p.refFE().phi ( j, fe_u.quadRule().quadPointCoor ( ig, 0 ), fe_u.quadRule().quadPointCoor ( ig, 1 ),
                                        fe_u.quadRule().quadPointCoor ( ig, 2 ) ) * fe_u.phiDer ( i, icoor, ig ) * fe_u.weightDet ( ig );
            mat ( i, j ) += coef * s;
        }
    }
}

void div ( const int icoor, Real coef, MatrixElemental& elmat,
           const CurrentFE& fe_u, const CurrentFE& fe_p,
           int iblock, int jblock )
/*
  \int q_i \frac{\partial v_j}{\partial x_icoor}
*/
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    UInt ig;
    UInt i, j;
    double s;
    for (i = 0; i < fe_u.nbFEDof(); i++)
    {
        for (j = 0; j < fe_p.nbFEDof(); j++)
        {
            s = 0;
            for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
            {
                s -= fe_u.phi (i, ig) * fe_p.phiDer (j, icoor, ig) * fe_u.weightDet ( ig );
            }
            mat ( i, j ) += coef * s;
        }
    }
}

void grad_div ( Real coef_grad, Real coef_div, MatrixElemental& elmat,
                const CurrentFE& fe_u, const CurrentFE& fe_p,
                int block_pres )
/*
  \int q_j \frac{\partial v_i}{\partial x_icoor}
*/
{
    double s;
    int iblock = block_pres - fe_u.nbLocalCoor();
    for ( UInt icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
    {
        MatrixElemental::matrix_view mat_grad = elmat.block ( iblock + icoor, block_pres );
        MatrixElemental::matrix_view mat_div = elmat.block ( block_pres , iblock + icoor );
        for ( UInt i = 0; i < fe_u.nbFEDof(); i++ )
        {
            for ( UInt j = 0; j < fe_p.nbFEDof(); j++ )
            {
                s = 0;
                for ( UInt ig = 0; ig < fe_u.nbQuadPt(); ig++ )
                {
                    s -= fe_p.phi ( j, ig ) * fe_u.phiDer ( i, icoor, ig ) * fe_u.weightDet ( ig );
                }
                mat_grad ( i, j ) += coef_grad * s;
                mat_div ( j, i ) += coef_div * s;
            }
        }
    }
}
//
void stab_stokes ( Real visc, Real coef_stab, MatrixElemental& elmat,
                   const CurrentFE& fe, int block_pres )
{
    MatrixElemental::matrix_view mat = elmat.block ( block_pres, block_pres );
    Real s, h = fe.diameter();
    Real fh2 = coef_stab * h * h / ( 2 * visc );
    for ( UInt i = 0; i < fe.nbFEDof(); i++ )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); j++ )
        {
            s = 0;
            for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
            {
                for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                {
                    s += fe.phiDer ( i, icoor, ig ) * fe.phiDer ( j, icoor, ig ) * fe.weightDet ( ig );
                }
            }
            mat ( i, j ) -= fh2 * s;
        }
    }
}

/*
 * Fixed by Umberto Villa,  Jan 2010
 */
void advection ( Real coef, VectorElemental& vel,
                 MatrixElemental& elmat, const CurrentFE& fe, int iblock, int jblock, int nb )
{
    Matrix mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );
    Real v_grad, s;
    Matrix v ( fe.nbQuadPt(), fe.nbLocalCoor() );

    //Evaluate the advective field at the quadrature nodes
    for ( UInt icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
    {
        VectorElemental::vector_view velicoor = vel.block ( icoor );
        for ( UInt iq = 0; iq < fe.nbQuadPt(); iq++ )
        {
            s = 0.;
            for ( UInt k = 0; k < fe.nbFEDof(); k++ )
            {
                s += velicoor ( k ) * fe.phi ( k, iq );    // velocity on the intgt point
            }
            v (iq, icoor) = s;
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
                for ( int icoor = 0; icoor < ( int ) fe.nbLocalCoor(); icoor++ )
                {
                    v_grad += v (iq, icoor) * fe.phiDer ( j, icoor, iq );
                }

                s += v_grad * fe.phi ( i, iq ) * fe.weightDet ( iq );
            }
            mat_tmp ( i, j ) = s * coef;
        }
    }

    // copy on the components
    for ( int icomp = 0; icomp < nb; icomp++ )
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block ( iblock + icomp, jblock + icomp );
        for ( UInt i = 0; i < fe.nbDiag(); i++ )
        {
            for ( UInt j = 0; j < fe.nbDiag(); j++ )
            {
                mat_icomp ( i, j ) += mat_tmp ( i, j );
            }
        }
    }
}

void grad ( const int icoor, const VectorElemental& vec_loc, MatrixElemental& elmat,
            const CurrentFE& fe1, const CurrentFE& fe2,
            int iblock, int jblock )
/*
  \int q_j \frac{\partial v_i}{\partial x_icoor}
*/
{
    //
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );

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
                    {
                        coef += vec_loc.vec() [ iloc + icoor * nbN1 ] * fe1.phi ( iloc, iq );
                    }

                    s += coef * fe2.phi ( i, iq ) * fe1.phiDer ( j, icoor, iq ) * fe1.weightDet ( iq );
                } // Loop on quadrature nodes

                mat ( i, j ) += s;
            } //Loop on j
        } // Loop on i
    } // if

}

//
// \! Gradient operator in the skew-symmetric form for NS Problems
//  A. Veneziani - December 2002
// \!
void grad_ss ( const int icoor, const VectorElemental& vec_loc, MatrixElemental& elmat,
               const CurrentFE& fe1, const CurrentFE& fe2,
               int iblock, int jblock )
/*
  \int vloc(icoor) \frac{\partial v_i}{\partial x_icoor} v_j + 1/2*\frac{\partial v_icoor}{\partial x_icoor} v_i v_j
*/
{
    //
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );

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
                        coef += vec_loc.vec() [ iloc + icoor * nbN1 ] * fe1.phi ( iloc, iq );
                        coef_div += vec_loc.vec() [ iloc + icoor * nbN1 ] * fe1.phiDer ( iloc, icoor, iq );
                    }

                    s += ( coef * fe1.phiDer ( j, icoor, iq ) + 0.5 * coef_div * fe1.phi ( j, iq ) ) * fe2.phi ( i, iq ) * fe1.weightDet ( iq );
                } // Loop on quadrature nodes

                mat ( i, j ) += s;
            } //Loop on j
        } // Loop on i
    } // if

}

// /!
// Gradient operator where the convective term is based on a local vector
// living on the basis given by fe3
// It is useful for advection diffusion problems driven by a NS problem
// !/
void grad ( const int icoor,
            const VectorElemental& vec_loc,
            MatrixElemental& elmat,
            const CurrentFE& fe1,
            const CurrentFE& fe2,
            const CurrentFE& fe3,
            int iblock, int jblock )
/*
  \int q_j \frac{\partial v_i}{\partial x_icoor}
*/
{
    //
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );

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
                    {
                        coef += vec_loc.vec() [iloc + icoor * nbN3] * fe3.phi (iloc, iq);
                    }

                    s += coef * fe2.phi (i, iq) * fe1.phiDer (j, icoor, iq) * fe1.weightDet (iq);
                } // Loop on quadrature nodes

                mat (i, j) += s;
            } //Loop on j
        } // Loop on i
    } // if

}


// Convective term with the velocity given in the quadrature nodes
void grad ( const int& icoor,
            const std::vector<Real>& localVector,
            MatrixElemental& elmat,
            const CurrentFE& currentFE1,
            const CurrentFE& currentFE2,
            const int& iblock,
            const int& jblock)
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );

    // This term concerns only the diagonal blocks (same components)
    if ( iblock == jblock )
    {

        for ( UInt iNode1 (0); iNode1 < currentFE1.nbFEDof(); iNode1++ )
        {
            for ( UInt jNode2 (0); jNode2 < currentFE2.nbFEDof(); jNode2++ )
            {
                Real sum (0.0);
                for ( UInt iq (0); iq < currentFE1.nbQuadPt(); iq++ )
                {
                    // Velocity in the quadrature node, component icoor
                    double coef (localVector[icoor * currentFE1.nbQuadPt() + iq]);

                    sum += coef
                           * currentFE2.phi (iNode1, iq)
                           * currentFE1.phiDer (jNode2, icoor, iq)
                           * currentFE1.weightDet (iq);
                } // Loop on quadrature nodes

                mat (iNode1, jNode2) += sum;
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
void quad(std::vector<Real> coef, MatrixElemental& elmat, VectorElemental& elvec,
const CurrentFE& fe,int iblock=0,int jblock=0)
{
MatrixElemental::matrix_view mat = elmat.block(iblock,jblock);
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
void source ( Real constant, VectorElemental& elvec, const CurrentFE& fe, int iblock )
{
    UInt i, ig;
    VectorElemental::vector_view vec = elvec.block ( iblock );
    Real s;
    for ( i = 0; i < fe.nbFEDof(); i++ )
    {
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
        {
            s += fe.phi ( i, ig ) * fe.weightDet ( ig );
        }
        vec ( i ) += constant * s;
    }
}

void source ( Real coef, VectorElemental& f, VectorElemental& elvec, const CurrentFE& fe,
              int fblock, int eblock )
/*
  compute \int f \phi_i
  where f is given on the dof of this element
  fblock is the block of the values read in f
  eblock is the block where the result is writen in the elvec
*/
{
    UInt i, ig;
    VectorElemental::vector_view vec = elvec.block ( eblock );
    VectorElemental::vector_view vecf = f.block ( fblock );
    Real f_ig;

    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        f_ig = 0.;
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {
            f_ig += vecf ( i ) * fe.phi ( i, ig );
        }
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {
            vec ( i ) += coef * f_ig * fe.phi ( i, ig ) * fe.weightDet ( ig );
        }
    }
}


void source_mass (const std::vector<Real>& constant, VectorElemental& elvec, const CurrentFE& currentFe, const int& iblock)
{
    VectorElemental::vector_view vec = elvec.block ( iblock );
    for (UInt iterNode (0); iterNode < currentFe.nbFEDof(); ++iterNode )
    {
        for ( UInt iterQuadNode (0); iterQuadNode < currentFe.nbQuadPt(); iterQuadNode++ )
        {
            vec (iterNode) += constant[iterQuadNode]
                              * currentFe.phi ( iterNode, iterQuadNode )
                              * currentFe.weightDet ( iterQuadNode );
        }
    }
}

void source_stiff (const std::vector<Real>& constant, VectorElemental& elvec, const CurrentFE& currentFe, const int& iblock)
{
    VectorElemental::vector_view vec = elvec.block ( iblock );
    const UInt nbQuadPt (currentFe.nbQuadPt() );
    for (UInt iterNode (0); iterNode < currentFe.nbFEDof(); ++iterNode )
    {
        for ( UInt iterQuadNode (0); iterQuadNode < nbQuadPt; iterQuadNode++ )
        {
            for (UInt iterGradComp (0); iterGradComp < currentFe.nbLocalCoor(); ++iterGradComp)
            {
                vec (iterNode) += constant[iterQuadNode + iterGradComp * nbQuadPt]
                                  * currentFe.phiDer ( iterNode, iterGradComp, iterQuadNode )
                                  * currentFe.weightDet ( iterQuadNode );
            }
        }
    }
}



void source_divuq (Real alpha, VectorElemental& uLoc,  VectorElemental& elvec, const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock  )
{
    UInt i, j, ic, iq;
    VectorElemental::vector_view vec = elvec.block ( iblock );
    Real s;

    for (i = 0; i < fe_p.nbFEDof(); i++)
    {
        s = 0;
        for (iq = 0; iq < fe_p.nbQuadPt(); ++iq )
            for (j = 0; j < fe_u.nbFEDof(); ++j)
                for (ic = 0; ic < fe_u.nbLocalCoor(); ++ic)
                {
                    s += uLoc[ic * fe_u.nbFEDof() + j] * fe_u.phiDer (j, ic, iq) * fe_p.phi (i, iq) * fe_p.weightDet ( iq );
                }

        vec ( i ) += s * alpha;
    }
}


void source_gradpv (Real alpha, VectorElemental& pLoc,  VectorElemental& elvec, const CurrentFE& fe_p, const CurrentFE& fe_u, int iblock )
{
    UInt i, j, iq;
    VectorElemental::vector_view vec = elvec.block ( iblock );
    Real s;

    for ( i = 0; i < fe_u.nbFEDof(); i++ )
    {
        s = 0;
        for (iq = 0; iq < fe_u.nbQuadPt(); ++iq )
            for (j = 0; j < fe_p.nbFEDof(); ++j)
            {
                s += pLoc[j] * fe_p.phiDer (j, iblock, iq) * fe_u.phi (i, iq) * fe_u.weightDet ( iq );
            }
        vec ( i ) += s * alpha;
    }
}




void source_fhn ( Real coef_f, Real coef_a, VectorElemental& u, VectorElemental& elvec, const CurrentFE& fe,
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
    VectorElemental::vector_view vec = elvec.block ( eblock );
    VectorElemental::vector_view vecu = u.block ( fblock );
    Real f_ig;

    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        f_ig = 0.;
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {
            f_ig += vecu ( i ) * fe.phi ( i, ig );
        }
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {
            vec ( i ) += coef_f * f_ig * ( 1 - f_ig ) * ( f_ig - coef_a ) * fe.phi ( i, ig ) * fe.weightDet ( ig );
        }
    }

}

//! \f$(beta\cdot\nabla u^k, v  )\f$
void source_advection ( const Real& coefficient, const VectorElemental& beta_loc, const VectorElemental& uk_loc,
                        VectorElemental& elvec, const CurrentFE& fe )
{
    boost::multi_array<Real, 2> guk (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()]);
    boost::multi_array<Real, 2> conv (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()]);
    std::vector<Real> beta (fe.nbLocalCoor() );

    Real s;

    UInt ig, icoor, jcoor, i;

    // loop on quadrature points
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // Interpolating beta
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {
            // Evaluate beta in the Gauss Point
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
            {
                s += fe.phi ( i, ig ) * beta_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
            beta[ icoor ] = s;
        }

        // Interpolation of grad u
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {
            // loop  on the derivative variable
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                // Evaluate the derivative in the gauss point
                s = 0.0;
                for ( i = 0; i < fe.nbFEDof(); i++ )
                {
                    //  grad u^k at a quadrature point
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
                }
                guk[ icoor ][ jcoor ] = s;
            }
        }

        // beta*(\grad u^k) at each quadrature point
        for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
        {
            s = 0.0;
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            {
                s += beta[ icoor ] * guk[ jcoor ][ icoor ];
            }
            conv[ ig ][ jcoor ] = s;
        }
    }

    //
    // Numerical integration
    //

    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
    {

        VectorElemental::vector_view vec = elvec.block ( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
            {
                s += conv[ ig ][ icoor ] * fe.phi ( i, ig ) * fe.wDetJacobian ( ig );
            }
            vec ( i ) += coefficient * s;
        }
    }
}

// coef * ( - \grad w^k :[I\div d - (\grad d)^T] u^k + convect^T[I\div d - (\grad d)^T] (\grad u^k)^T , v  ) for Newton FSI
//
// Remark: convect = u^n-w^k relative vel.
//
void source_mass1 ( Real coef, const VectorElemental& uk_loc, const VectorElemental& wk_loc, const VectorElemental& convect_loc,
                    const VectorElemental& d_loc, VectorElemental& elvec, const CurrentFE& fe )
{

    boost::multi_array<Real, 2> A (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()]);
    boost::multi_array<Real, 2> B (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()]);
    boost::multi_array<Real, 2> uk (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()]);
    boost::multi_array<Real, 3> guk (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()][fe.nbLocalCoor()]);
    std::vector<Real> dw (fe.nbLocalCoor() );
    std::vector<Real> aux (fe.nbQuadPt() );
    std::vector<Real> convect (fe.nbLocalCoor() );
    boost::multi_array<Real, 2> convect_A (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()]);

    Real s, sA, sB, sG;


    UInt icoor, jcoor, ig, i;
    // loop on quadrature points
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordindates
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // each compontent of uk at each quadrature points
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
            {
                s += fe.phi ( i, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
            uk[ ig ][ icoor ] = s;//uk_x(pt_ig), uk_y(pt_ig), uk_z(pt_ig)

            // each compontent of convect at this quadrature point
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
            {
                s += fe.phi ( i, ig ) * convect_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
            convect[ icoor ] = s;


            // loop  on space coordindates
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                sB = 0.0;
                sA = 0.0;
                sG = 0.0;
                for ( i = 0; i < fe.nbFEDof(); i++ )
                {
                    sG += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at each quadrature point
                    sB -= fe.phiDer ( i, jcoor, ig ) * wk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad (- w^k) at this quadrature point
                    sA -= fe.phiDer ( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe.nbFEDof() ]; //  - (\grad d) ^T at this quadrature point
                }
                guk[ ig ][ icoor ][ jcoor ] = sG; // \grad u^k at each quadrature point
                B[ icoor ][ jcoor ] = sB; // \grad (convect) at this quadrature point
                A[ icoor ][ jcoor ] = sA; // -(\grad d) ^T at this quadrature point
            }
        }

        s = 0.0;
        for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
        {
            s -= A[ jcoor ][ jcoor ];    // \div d at this quadrature point ( - trace( A ) )
        }

        for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
        {
            A[ jcoor ][ jcoor ] += s;    // I\div d - (\grad d)^T at this quadrature point (A+I(-tr(A)))
        }

        s = 0;
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s += B[ icoor ][ jcoor ] * A[ icoor ][ jcoor ];    // \grad (-w^k):[I\div d - (\grad d)^T] at each quadrature point
            }
        aux[ ig ] = s;

        s = 0;
        for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
        {
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            {
                s += convect[ icoor ] * A[ icoor ][ jcoor ];    // convect^T [I\div d - (\grad d)^T]
            }
            convect_A[ ig ][ jcoor ] = s;
        }
        //         std::cout<<"aux "<<aux[ig]<<std::endl;

        //     for(icoor=0; icoor<fe.nbLocalCoor(); ++icoor)
        //         for(jcoor=0; jcoor<fe.nbLocalCoor(); ++jcoor)
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
    for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
    {

        // the block iccor of the elementary vector
        VectorElemental::vector_view vec = elvec.block ( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
            {

                // \grad ( - w^k ):[I\div d - (\grad d)^T] \phi_i
                s += aux[ ig ] * uk[ ig ][ icoor ] * fe.phi ( i, ig ) * fe.weightDet ( ig );

                // convect^T [I\div d - (\grad d)^T] (\grad u^k)^T \phi_i
                for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                {
                    s += convect_A[ ig ][ jcoor ] * guk[ ig ][ icoor ][ jcoor ] * fe.phi ( i, ig ) * fe.weightDet ( ig );
                }
            }
            vec ( i ) += coef * s;
        }
    }
}



//
// coef * ( \grad u^k dw, v  ) for Newton FSI
//
//
void source_mass2 ( Real coef, const VectorElemental& uk_loc, const VectorElemental& dw_loc,
                    VectorElemental& elvec, const CurrentFE& fe )
{

    boost::multi_array<Real, 2> guk (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()]);
    std::vector<Real> dw (fe.nbLocalCoor() );
    boost::multi_array<Real, 2> aux (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()]);

    Real s;

    UInt ig, icoor, jcoor, i;

    // loop on quadrature points
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // each compontent (icoor) of dw at this quadrature point
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
            {
                s += fe.phi ( i, ig ) * dw_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
            dw[ icoor ] = s;

            // loop  on space coordinates
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s = 0.0;
                for ( i = 0; i < fe.nbFEDof(); i++ )
                {
                    s += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];    //  \grad u^k at a quadrature point
                }
                guk[ icoor ][ jcoor ] = s;
            }
        }

        // (\grad u^k)dw at each quadrature point
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {
            s = 0.0;
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s += guk[ icoor ][ jcoor ] * dw[ jcoor ];
            }
            aux[ ig ][ icoor ] = s;
        }
    }

    //
    // Numerical integration
    //

    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
    {

        VectorElemental::vector_view vec = elvec.block ( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
            {
                s += aux[ ig ][ icoor ] * fe.phi ( i, ig ) * fe.weightDet ( ig );
            }
            vec ( i ) += coef * s;
        }
    }
}



//
// coef * ( \grad u^n :[2 I \div d - (\grad d)^T]  u^k , v  ) for Newton FSI
//
//
void source_mass3 ( Real coef, const VectorElemental& un_loc, const VectorElemental& uk_loc, const VectorElemental& d_loc,
                    VectorElemental& elvec, const CurrentFE& fe )
{

    boost::multi_array<Real, 2> B (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()]);
    boost::multi_array<Real, 2> A (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()]);
    boost::multi_array<Real, 2> uk (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()]);
    std::vector<Real> aux (fe.nbQuadPt() );

    Real s, sA, sB;


    UInt icoor, jcoor, ig, i;
    // loop on quadrature points
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {

        // loop on space coordindates
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // each compontent of uk at each quadrature points
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
            {
                s += fe.phi ( i, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
            uk[ ig ][ icoor ] = s;


            // loop  on space coordindates
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                sB = 0.0;
                sA = 0.0;
                for ( i = 0; i < fe.nbFEDof(); i++ )
                {
                    sB += fe.phiDer ( i, jcoor, ig ) * un_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^n at this quadrature point
                    sA -= fe.phiDer ( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe.nbFEDof() ]; //  - (\grad d) ^T at this quadrature point
                }
                B[ icoor ][ jcoor ] = sB; // \grad u^n at this quadrature point
                A[ icoor ][ jcoor ] = sA; // -(\grad d) ^T at this quadrature point
            }
        }

        s = 0.0;
        for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
        {
            s -= A[ jcoor ][ jcoor ];    // \div d at this quadrature point ( - trace( A ) )
        }

        for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
        {
            A[ jcoor ][ jcoor ] += 2 * s;    // 2 * I\div d - (\grad d)^T at this quadrature point
        }

        s = 0;
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s += B[ icoor ][ jcoor ] * A[ icoor ][ jcoor ];    // \grad u^n:[2 * I\div d - (\grad d)^T] at each quadrature point
            }
        aux[ ig ] = s;
    }

    // At this point we have:
    //    v u^k at each quadrature point: uk
    //    v  \grad u^n:[ 2 * I\div d - (\grad d)^T]: aux


    //
    // Numerical integration
    //

    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
    {

        // the block iccor of the elementary vector
        VectorElemental::vector_view vec = elvec.block ( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {
            // loop on quadrature points
            s = 0;
            for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
                // \grad u^n:[2 * I\div d - (\grad d)^T] u^k \phi_i
            {
                s += aux[ ig ] * uk[ ig ][ icoor ] * fe.phi ( i, ig ) * fe.weightDet ( ig );
            }
            vec ( i ) += coef * s;
        }
    }
}






//
// coef * ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  ) for Newton FSI
//
void source_stress ( Real coef, Real mu, const VectorElemental& uk_loc, const VectorElemental& pk_loc,
                     const VectorElemental& d_loc, VectorElemental& elvec, const CurrentFE& fe_u,
                     const CurrentFE& fe_p )
{

    boost::multi_array<Real, 3> B (
        boost::extents[fe_u.nbQuadPt()][fe_u.nbLocalCoor()][fe_u.nbLocalCoor()]);
    boost::multi_array<Real, 2> A (
        boost::extents[fe_u.nbLocalCoor()][fe_u.nbLocalCoor()]);
    boost::multi_array<Real, 2> guk (
        boost::extents[fe_u.nbLocalCoor()][fe_u.nbLocalCoor()]);
    boost::multi_array<Real, 2> sigma (
        boost::extents[fe_u.nbLocalCoor()][fe_u.nbLocalCoor()]);

    Real s, sA, sG, pk;

    UInt icoor, jcoor, kcoor, ig, i;

    // loop on quadrature points
    for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordindates
            for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
            {
                sA = 0.0;
                sG = 0.0;
                for ( i = 0; i < fe_u.nbFEDof(); i++ )
                {
                    sG += fe_u.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbFEDof() ]; //  \grad u^k at this quadrature point
                    sA -= fe_u.phiDer ( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe_u.nbFEDof() ]; //  - (\grad d) ^T at this quadrature point
                }
                guk[ icoor ][ jcoor ] = sG;
                A[ icoor ][ jcoor ] = sA;
            }
        }

        s = 0.0;
        for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
        {
            s -= A[ jcoor ][ jcoor ];    // \div d at a quadrature point ( - trace( A ) )
        }

        //         std::cout<<"div = "<< s <<std::endl;

        for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
        {
            A[ jcoor ][ jcoor ] += s;    // I\div d  - (\grad d)^T
        }

        pk = 0.0;
        for ( i = 0; i < fe_p.nbFEDof(); i++ )
        {
            pk += fe_p.phi ( i, ig ) * pk_loc.vec() [ i ];    // p^k at this quadrature point
        }


        // sigma = [-p^k I + 2*mu e(u^k)] a quadrature point
        for ( icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
        {
            for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
            {
                sigma[ icoor ][ jcoor ] = mu * ( guk[ icoor ][ jcoor ] + guk[ jcoor ][ icoor ] );
            }
            sigma[ icoor ][ icoor ] -= pk;
        }

        // [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] at each quadrature point
        for ( icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
            for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
            {
                s = 0;
                for ( kcoor = 0; kcoor < fe_u.nbLocalCoor(); kcoor++ )
                {
                    s += sigma[ icoor ][ kcoor ] * A[ kcoor ][ jcoor ];
                }
                B[ ig ][ icoor ][ jcoor ] = s;
            }
    }

    //     for(icoor=0; icoor<fe_u.nbLocalCoor(); ++icoor)
    //         for(jcoor=0; jcoor<fe_u.nbLocalCoor(); ++jcoor)
    //             {
    //                 std::cout<<" Atrue ["<<icoor<<"]["<<jcoor<<"] = "<<A[icoor][jcoor]<<std::endl;
    //             }

    //     for(icoor=0; icoor<fe_u.nbLocalCoor(); ++icoor)
    //         for(jcoor=0; jcoor<fe_u.nbLocalCoor(); ++jcoor)
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
    for ( icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
    {

        VectorElemental::vector_view vec = elvec.block ( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0; i < fe_u.nbFEDof(); i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )

                for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
                {
                    s += B[ ig ][ icoor ][ jcoor ] * fe_u.phiDer ( i, jcoor, ig ) * fe_u.weightDet ( ig );
                }
            vec ( i ) += coef * s;
        }
    }
}



//
// + \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )
//
void source_stress2 ( Real coef, const VectorElemental& uk_loc, const VectorElemental& d_loc, VectorElemental& elvec, const CurrentFE& fe_u )
{

    boost::multi_array<Real, 3> A (
        boost::extents[fe_u.nbQuadPt()][fe_u.nbLocalCoor()][fe_u.nbLocalCoor()]);
    boost::multi_array<Real, 2> guk (
        boost::extents[fe_u.nbLocalCoor()][fe_u.nbLocalCoor()]);
    boost::multi_array<Real, 2> gd (
        boost::extents[fe_u.nbLocalCoor()][fe_u.nbLocalCoor()]);

    Real su, sd, s;

    UInt icoor, jcoor, kcoor, ig, i;

    // loop on quadrature points
    for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordindates
            for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
            {
                su = 0.0;
                sd = 0.0;
                for ( i = 0; i < fe_u.nbFEDof(); i++ )
                {
                    su += fe_u.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbFEDof() ]; //  \grad u^k at this quadrature point
                    sd += fe_u.phiDer ( i, jcoor, ig ) * d_loc.vec() [ i + icoor * fe_u.nbFEDof() ]; //  \grad d at this quadrature point
                }
                guk[ icoor ][ jcoor ] = su;
                gd[ icoor ][ jcoor ] = sd;
            }
        }


        for ( icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
        {
            for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
            {
                s = 0;
                for ( kcoor = 0; kcoor < fe_u.nbLocalCoor(); kcoor++ )
                {
                    s += guk[ icoor ][ kcoor ] * gd[ kcoor ][ jcoor ] + gd[ kcoor ][ icoor ] * guk[ jcoor ][ kcoor ];
                }
                A[ ig ][ icoor ][ jcoor ] = s;
            }
        }

        //     for(icoor=0; icoor<fe_u.nbLocalCoor(); ++icoor)
        //         for(jcoor=0; jcoor<fe_u.nbLocalCoor(); ++jcoor)
        //             {
        //                 std::cout<<" gdtrue ["<<icoor<<"]["<<jcoor<<"] = "<<gd[icoor][jcoor]<<std::endl;
        //             }
        //     for(icoor=0; icoor<fe_u.nbLocalCoor(); ++icoor)
        //         for(jcoor=0; jcoor<fe_u.nbLocalCoor(); ++jcoor)
        //             {
        //                 std::cout<<" Atrue ["<<icoor<<"]["<<jcoor<<"] = "<<A[ig][icoor][jcoor]<<std::endl;
        //             }
    }

    //
    // Numerical integration
    //
    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
    {

        VectorElemental::vector_view vec = elvec.block ( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0; i < fe_u.nbFEDof(); i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )

                for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
                {
                    s += fe_u.phiDer ( i, jcoor, ig ) * A[ ig ][ icoor ][ jcoor ] * fe_u.weightDet ( ig );
                }
            vec ( i ) += coef * s;
        }
    }
}





//
// coef * (  (\grad u^k):[I\div d - (\grad d)^T] , q  ) for Newton FSI
//
void source_press ( Real coef, const VectorElemental& uk_loc, const VectorElemental& d_loc, VectorElemental& elvec,
                    const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock )
{
    boost::multi_array<Real, 2> A (
        boost::extents[fe_u.nbLocalCoor()][fe_u.nbLocalCoor()]);
    boost::multi_array<Real, 2> guk (
        boost::extents[fe_u.nbLocalCoor()][fe_u.nbLocalCoor()]);
    std::vector<Real> aux (fe_u.nbQuadPt() );

    VectorElemental::vector_view vec = elvec.block ( iblock );

    Real s, sA, sG;
    UInt icoor, jcoor, ig, i;


    // loop on quadrature points
    for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
        {

            // loop  on space coordinates
            for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
            {
                sA = 0.0;
                sG = 0.0;
                for ( i = 0; i < fe_u.nbFEDof(); i++ )
                {
                    sG += fe_u.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbFEDof() ]; //  \grad u^k at a quadrature point
                    sA -= fe_u.phiDer ( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe_u.nbFEDof() ]; //  - (\grad d) ^T at a quadrature point
                }
                guk[ icoor ][ jcoor ] = sG;
                A[ icoor ][ jcoor ] = sA;
            }
        }

        s = 0.0;
        for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
        {
            s -= A[ jcoor ][ jcoor ];    // \div d at this quadrature point ( - trace( A ) )
        }

        for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
        {
            A[ jcoor ][ jcoor ] += s;    // I\div d - (\grad d)^T at this quadrature point
        }

        s = 0;
        for ( icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
            for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
            {
                s += guk[ icoor ][ jcoor ] * A[ icoor ][ jcoor ];    // \grad u^k : [I\div d - (\grad d)^T] at each quadrature point
            }
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
            s += aux[ ig ] * fe_p.phi ( i, ig ) * fe_u.weightDet ( ig );
        }
        vec [ i ] += coef * s;
    }
}


//
// coef * ( [I\div d - (\grad d)^T - \grad d] \grap p, \grad q  ) for Newton FSI
//
void source_press2 ( Real coef, const VectorElemental& p_loc, const VectorElemental& d_loc, VectorElemental& elvec,
                     const CurrentFE& fe, int iblock )
{
    boost::multi_array<Real, 2> A (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()]);
    boost::multi_array<Real, 2> B (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()]);
    boost::multi_array<Real, 2> aux (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()]);
    std::vector<Real> gpk (fe.nbLocalCoor() );

    VectorElemental::vector_view vec = elvec.block ( iblock );

    Real s, sA, sG;
    UInt icoor, jcoor, ig, i;


    // loop on quadrature points
    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {


        // loop on space coordinates
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            sG = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
            {
                sG += fe.phiDer ( i, icoor, ig ) * p_loc.vec() [ i ];    //  \grad p^k at a quadrature point
            }
            gpk[ icoor ] = sG;

            // loop  on space coordinates
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                sA = 0.0;
                for ( i = 0; i < fe.nbFEDof(); i++ )
                {
                    sA -= fe.phiDer ( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe.nbFEDof() ];    //  - (\grad d) ^T at a quadrature point
                }
                A[ icoor ][ jcoor ] = sA;
                B[ jcoor ][ icoor ] = sA;
            }
        }

        s = 0.0;
        for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
        {
            s -= A[ jcoor ][ jcoor ];    // \div d at this quadrature point ( - trace( A ) )
        }

        for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
        {
            A[ jcoor ][ jcoor ] += s;    // I\div d - (\grad d)^T at this quadrature point
        }

        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                A[ icoor ][ jcoor ] += B[ icoor ][ jcoor ];    // I\div d - (\grad d)^T - \grad d at this quadrature point
            }

        s = 0;
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                s += A[ icoor ][ jcoor ] * gpk[jcoor];    // [I\div d - (\grad d)^T -\grad d] \grad p^k at each quadrature point
            }
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
            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            {
                s += aux[ ig ][icoor] * fe.phiDer ( i, icoor, ig ) * fe.weightDet ( ig );
            }
        vec [ i ] += coef * s;
    }
}





//
//Cholesky decomposition
void choldc ( KNM<Real>& a, KN<Real>& p )
{
    int i, j, k;
    Real sum;

    int n = a.N();
    for ( i = 0; i < n; i++ )
    {
        for ( j = i; j < n; j++ )
        {
            for ( sum = a ( i, j ), k = i - 1; k >= 0; k-- )
            {
                sum -= a ( i, k ) * a ( j, k );
            }
            if ( i == j )
            {
                p ( i ) = std::sqrt ( sum );
            }
            else
            {
                a ( j, i ) = sum / p ( i );
            }
        }
    }
}
//
//Cholesky solution
void cholsl ( KNM<Real>& a, KN<Real>& p, KN<Real>& b, KN<Real>& x )
{
    int i, k;
    Real sum;

    int n = a.N();
    for ( i = 0; i < n; i++ )
    {
        for ( sum = b ( i ), k = i - 1; k >= 0; k-- )
        {
            sum -= a ( i, k ) * x ( k );
        }
        x ( i ) = sum / p ( i );
    }
    for ( i = n - 1; i >= 0; i-- )
    {
        for ( sum = x ( i ), k = i + 1; k < n; k++ )
        {
            sum -= a ( k, i ) * x ( k );
        }
        x ( i ) = sum / p ( i );
    }
}

//
// coef * (  (\grad u^k):[I\div d - (\grad d)^T] , q  ) for Newton FSI
//
void source_press ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat,
                    const CurrentFE& fe_u, const CurrentFE& fe_p, ID mmDim , int iblock )
{
    boost::multi_array<Real, 4> A (
        boost::extents[fe_u.nbLocalCoor()][fe_u.nbLocalCoor()][fe_u.nbFEDof()][fe_u.nbLocalCoor()]);
    boost::multi_array<Real, 2> guk (
        boost::extents[fe_u.nbLocalCoor()][fe_u.nbLocalCoor()]);
    boost::multi_array<Real, 3> aux (
        boost::extents[fe_u.nbQuadPt()][fe_u.nbFEDof()][fe_u.nbLocalCoor()]);

    Real l, sA, sG;
    UInt icoor, jcoor, ig;


    // loop on quadrature points
    for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
    {


        ////INIT
        for (UInt p = 0; p < fe_u.nbLocalCoor(); ++p)
        {
            for (UInt q = 0; q < fe_u.nbLocalCoor(); ++q)
            {
                for (UInt d = 0; d < fe_u.nbFEDof(); ++d)
                {
                    for (UInt e = 0; e < fe_u.nbLocalCoor(); ++e)
                    {
                        A[p][q][d][e] = 0.;
                    }
                }
                //guk[p][q]=0.;
            }
            //for(int h=0; h<fe_u.nbFEDof(); ++h)
            //aux[ig][h][p]=0.;
        }
        ////END INIT

        // loop on space coordinates
        for ( icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
        {

            //for ( i = 0;i < fe_u.nbFEDof();i++ )
            //for ( short k = 0;k < fe_u.nbLocalCoor();k++ )
            //sA[i][k] = 0.0;
            // loop  on space coordinates
            for ( UInt kcoor = 0; kcoor < fe_u.nbLocalCoor(); kcoor++ )
                for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
                {
                    sG = 0.0;
                    for ( UInt i = 0; i < fe_u.nbFEDof(); i++ )
                    {
                        sG += fe_u.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbFEDof() ]; //  \grad u^k at a quadrature point
                        sA = -fe_u.phiDer ( i, icoor, ig ) /**fe_u.phi( i, ig )*/ ;/** d_loc.vec() [ i + jcoor * fe_u.nbFEDof() ];*/ //  - (\grad d) ^T at a quadrature point
                        A[ icoor ][ jcoor ][i][jcoor] = sA; // -\delta_{jcoor kcoor} \partial_{icoor}
                    }
                    guk[ icoor ][ jcoor ] = sG;
                }
        }





        std::vector<Real> z (fe_u.nbLocalCoor() );
        for ( UInt i = 0; i < fe_u.nbFEDof(); i++ )
        {
            //                for ( int kcoor = 0;kcoor < fe_u.nbLocalCoor();kcoor++ )
            for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
            {
                //if(icoor==jcoor)
                z[ jcoor ] = A[ jcoor ][ jcoor ][i][jcoor];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
            }

            for ( jcoor = 0; jcoor < fe_p.nbLocalCoor(); jcoor++ )
            {
                for (UInt kcoor = 0; kcoor < fe_p.nbLocalCoor(); kcoor++ )
                {
                    //if(icoor==jcoor)
                    A[ jcoor ][ jcoor ][i][kcoor] -= z[ kcoor ];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
                }
            }



            //                         for ( jcoor = 0;jcoor < fe_u.nbLocalCoor();jcoor++ )
            //                             {
            //                             //if(kcoor==jcoor)
            //                                     for ( icoor = 0;icoor < fe_u.nbLocalCoor();icoor++ )
            //                                         {
            //                                             //    if(icoor==jcoor)
            //                                             A[ icoor ][ jcoor ][i][kcoor] -= A[ kcoor ][ icoor ][i][jcoor];  // \delta_{jcoor icoor}\partial_{kcoor}
            //                                             //else
            //                                             //f[icoor][jcoor]=0.;
            //                                         }
            //                             }
            //                    }
            //                for ( UInt kcoor = 0;kcoor < fe_u.nbLocalCoor();kcoor++ )
            //                    {
            //                         for ( jcoor = 0;jcoor < fe_u.nbLocalCoor();jcoor++ )
            //                         {
            //                             for ( icoor = 0;icoor < fe_u.nbLocalCoor();icoor++ )
            //                                 {
            //                                     //for ( jcoor = 0;jcoor < fe_u.nbLocalCoor();jcoor++ )
            //                                     if(icoor==jcoor)
            //                                         A[ icoor ][ jcoor ][i][kcoor] += f[icoor][jcoor];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
            //                 //                for ( i = 0;i < fe_u.nbFEDof();i++ )
            //                 //                    for(jcoor=0;jcoor<fe_u.nbLocalCoor();++jcoor)
            //                 //                        s[i][jcoor] = 0.0;
            //                               //for ( jcoor = 0;jcoor < fe_u.nbLocalCoor();jcoor++ )
            //                               //{
            //                                 }
            //                         }
            //                    }
            //                for ( UInt kcoor = 0;kcoor < fe_u.nbLocalCoor();kcoor++ )
            //                    {
            for ( UInt kcoor = 0; kcoor < fe_u.nbLocalCoor(); kcoor++ )
            {
                l = 0;
                //for ( kcoor = 0;kcoor < fe_u.nbLocalCoor();kcoor++ )
                for ( icoor = 0; icoor < fe_u.nbLocalCoor(); icoor++ )
                    for ( jcoor = 0; jcoor < fe_u.nbLocalCoor(); jcoor++ )
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
    for ( UInt kcoor = 0; kcoor < fe_u.nbLocalCoor(); kcoor++ )
    {
        //            for ( short l = 0;i < fe_u.nbFEDof();i++ )
        //for ( jcoor = 0;jcoor < fe_u.nbLocalCoor();jcoor++ )
        double l = 0.;

        MatrixElemental::matrix_view mat = elmat.block ( iblock, kcoor );
        for ( UInt j = 0; j < mmDim; j++ )
            for ( UInt i = 0; i < fe_p.nbFEDof(); i++ )
            {
                mat (i, j) = 0.;
            }

        // Loop on nodes, i.e. loop on elementary vector components
        for ( UInt j = 0; j < mmDim; j++ )
            for (UInt i = 0; i < fe_p.nbFEDof(); i++ )
            {
                l = 0.;
                // loop on quadrature points
                for ( ig = 0; ig < fe_u.nbQuadPt(); ig++ )
                {
                    //            std::cout << ig << " " << fe_p.phi(i, ig) << std::endl;
                    l += aux[ ig ][j][kcoor] * fe_p.phi ( i, ig ) * fe_u.weightDet ( ig );
                }
                mat ( i , j) += coef * l;
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
void shape_terms (
    //const VectorElemental& d_loc,
    Real rho,
    Real mu,
    const VectorElemental& un_loc,
    const VectorElemental& uk_loc,
    const VectorElemental& wk_loc,
    const VectorElemental& convect_loc,
    const VectorElemental& pk_loc,
    MatrixElemental& elmat,
    const CurrentFE& fe,
    const CurrentFE& fe_p,
    ID /*mmDim*/,
    MatrixElemental& /*elmatP*/,
    int /*iblock*/,
    bool wImplicit,
    Real alpha,
    boost::shared_ptr<MatrixElemental> elmat_convect
)
{
    // I div d - (grad d)^T
    boost::multi_array<Real, 4> eta (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbFEDof()][fe.nbLocalCoor()]);
    // I div d - (grad d)^T
    boost::multi_array<Real, 4> etaMass3 (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbFEDof()][fe.nbLocalCoor()]);
    // u^k
    boost::multi_array<Real, 2> uk (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()]);
    // grad u^k
    boost::multi_array<Real, 3> guk (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()][fe.nbLocalCoor()]);
    // div u^k
    std::vector<Real> duk (fe.nbQuadPt() );
    // grad u^k
    boost::multi_array<Real, 3> gun (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()][fe.nbLocalCoor()]);
    // convect
    std::vector<Real> convect (fe.nbLocalCoor() );
    // convect^T [ I div d - (grad d)^T ]
    boost::multi_array<Real, 4> convect_eta (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()][fe.nbFEDof()][fe.nbLocalCoor()]);
    // - p^k I + 2 mu e(u^k)
    boost::multi_array<Real, 2> sigma (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()]);
    // ( - p^k I + 2 mu e(u^k) )( I div d - (grad d)^T )
    boost::multi_array<Real, 5> B (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbFEDof()][fe.nbLocalCoor()]);
    // gd
    boost::multi_array<Real, 4> gd (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbFEDof()][fe.nbLocalCoor()]);
    // grad u^k grad d + (grad d)^T (grad u^k)^T
    boost::multi_array<Real, 5> A2 (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()][fe.nbLocalCoor()][fe.nbFEDof()][fe.nbLocalCoor()]);

    boost::multi_array<Real, 4> aux (
        boost::extents[fe.nbQuadPt()][fe.nbLocalCoor()][fe.nbFEDof()][fe.nbLocalCoor()]);
    boost::multi_array<Real, 2> BMass (
        boost::extents[fe.nbLocalCoor()][fe.nbLocalCoor()]);
    boost::multi_array<Real, 3> auxMass (
        boost::extents[fe.nbQuadPt()][fe.nbFEDof()][fe.nbLocalCoor()]);
    boost::multi_array<Real, 3> auxMass3 (
        boost::extents[fe.nbQuadPt()][fe.nbFEDof()][fe.nbLocalCoor()]);

    Real s,  sA, sB, sGk, sGn, pk;

    UInt icoor, jcoor, ig, kcoor, i, j;

    // loop on quadrature points

    for ( ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        ////INIT
        for (UInt p = 0; p < fe.nbLocalCoor(); ++p)
        {
            uk[ig][p] = 0.;
            convect[p] = 0.;
            for (UInt q = 0; q < fe.nbLocalCoor(); ++q)
            {
                guk[ig][p][q] = 0.;
                gun[ig][p][q] = 0.;
                sigma[p][q] = 0.;
                BMass[p][q] = 0.;
                for (UInt d = 0; d < fe.nbFEDof(); ++d)
                {
                    auxMass[ ig ][d][q] = 0.;
                    auxMass3[ ig ][d][q] = 0.;
                    aux[ig][p][d][q] = 0.;
                    convect_eta[ig][p][d][q] = 0.;
                    for (UInt e = 0; e < fe.nbLocalCoor(); ++e)
                    {
                        gd[p][q][d][e] = 0.;
                        eta[p][q][d][e] = 0.;
                        etaMass3[p][q][d][e] = 0.;
                        A2[ig][p][q][d][e] = 0.;
                        B[ig][p][q][d][e] = 0.;
                    }
                }
            }
        }
        ////////END INIT

        // loop on space coordindates
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {

            // each compontent of uk at each quadrature points
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
            {
                s += fe.phi ( i, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
            uk[ ig ][ icoor ] = s;//uk_x(pt_ig), uk_y(pt_ig), uk_z(pt_ig)

            // each compontent of convect at this quadrature point
            s = 0.0;
            for ( i = 0; i < fe.nbFEDof(); i++ )
            {
                s += fe.phi ( i, ig ) * convect_loc.vec() [ i + icoor * fe.nbFEDof() ];
            }
            convect[ icoor ] = s;


            // loop  on space coordindates
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                sB = 0.0;
                sGk = 0.0;
                sGn = 0.0;
                for ( i = 0; i < fe.nbFEDof(); i++ )
                {
                    gd[ icoor ][ jcoor ][i][jcoor] = fe.phiDer ( i, jcoor, ig ) /** d_loc.vec() [ i + jcoor * fe.nbFEDof() ]*/;

                    sGk += fe.phiDer ( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at each quadrature point
                    sGn += fe.phiDer ( i, jcoor, ig ) * un_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad u^k at each quadrature point
                    sB -= fe.phiDer ( i, jcoor, ig ) * wk_loc.vec() [ i + icoor * fe.nbFEDof() ]; //  \grad (- w^k) at this quadrature point
                    sA = -fe.phiDer ( i, icoor, ig ) /** d_loc.vec() [ i + jcoor * fe.nbFEDof() ]*/; //  - (\grad d) ^T at this quadrature point
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
        {
            pk += fe_p.phi ( i, ig ) * pk_loc.vec() [ i ];    // p^k at this quadrature point
        }

        // sigma = [-p^k I + 2*mu e(u^k)] a quadrature point
        for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
        {
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                sigma[ icoor ][ jcoor ] = mu * ( guk[ig][ icoor ][ jcoor ] + guk[ig][ jcoor ][ icoor ] );
            }
            sigma[ icoor ][ icoor ] -= pk;
        }

        //!building the tensor \f$\eta = [I\nabla\cdot d - (\nabla d)^T]\f$
        for ( i = 0; i < fe.nbFEDof(); i++ )
        {

            std::vector<Real> z (fe.nbLocalCoor() );
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                //if(icoor==jcoor)
                z[ jcoor ] = eta[ jcoor ][ jcoor ][i][jcoor];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
            }
            for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
            {
                for ( kcoor = 0; kcoor < fe.nbLocalCoor(); kcoor++ )
                {
                    eta[ jcoor ][ jcoor ][i][kcoor] -= z[kcoor];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
                }
            }

            //!source_mass1

            for ( kcoor = 0; kcoor < fe.nbLocalCoor(); kcoor++ )
            {
                s = 0;
                for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                    for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                    {
                        s += BMass[ icoor ][ jcoor ] * eta[ icoor ][ jcoor ][i][kcoor];    // \grad (-w^k):[I\div d - (\grad d)^T] at each quadrature point
                    }
                auxMass[ ig ][i][kcoor] = s;

                for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                {
                    s = 0.;
                    for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                    {
                        s += convect[ icoor ] * eta[ icoor ][ jcoor ][i][kcoor]; // convect^T [I\div d - (\grad d)^T]
                    }
                    convect_eta[ ig ][ jcoor ][i][kcoor] = s;
                }
            }
            //! source_stress

            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
            {
                for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                {
                    for ( kcoor = 0; kcoor < fe.nbLocalCoor(); kcoor++ )
                    {
                        s = 0.;
                        for (UInt zcoor = 0; zcoor < fe.nbLocalCoor(); zcoor++ )
                        {
                            s += sigma[ icoor ][ zcoor ] * eta[ zcoor ][ jcoor ][i][kcoor];
                        }
                        B[ ig ][ icoor ][ jcoor ][i][kcoor] = s;
                    }
                }
            }


            //! source_stress2
            for ( kcoor = 0; kcoor < fe.nbLocalCoor(); kcoor++ )
            {
                for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                {
                    for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                    {
                        s = 0.;
                        for ( UInt zcoor = 0; zcoor < fe.nbLocalCoor(); zcoor++ )
                            // \grad u^k \grad d + [\grad d]^T[\grad u^k]^T  at each quadrature point
                        {
                            s += guk[ig][ icoor ][ zcoor ] * gd[ zcoor ][ jcoor ][i][kcoor] + gd[ zcoor ][ icoor ][i][kcoor] * guk[ig][ jcoor ][ zcoor ];
                        }
                        A2[ ig ][ icoor ][ jcoor ][i][kcoor] = s;
                    }
                }
            }

            if (wImplicit) //source_mass2 and convective term derivative
            {
                for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                {
                    //s = 0.0;
                    for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                    {
                        aux[ ig ][ icoor ][i][jcoor] = guk[ig][ icoor ][ jcoor ]  * fe.phi ( i, ig ) /**alpha*/;
                    }

                }
                for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                {
                    duk[ig] += guk[ig][ jcoor ][ jcoor ];  //div uk
                }
            }


            //source_mass3
            // coef * ( \grad u^n :[2 I \div d - (\grad d)^T]  u^k , v  ) for Newton FSI
            //
            //

            //!building the tensor \f$\eta = [I\nabla\cdot d - (\nabla d)^T]\f$
            //                for ( int kcoor = 0;kcoor < fe.nbLocalCoor();kcoor++ )

            for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                {
                    for (kcoor = 0; kcoor < fe.nbLocalCoor(); kcoor++ )
                    {
                        //if(icoor==jcoor)
                        etaMass3[ icoor ][ jcoor ][i][kcoor] -= 2 * z[kcoor]; //! -2*\delta_{jcoor, kcoor} \partial_{icoor} + \delta_{jcoor, icoor}\partial_{kcoor}
                    }
                }

            for (kcoor = 0; kcoor < fe.nbLocalCoor(); ++kcoor)
            {
                s = 0;
                for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
                    for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                    {
                        s += gun[ig][ icoor ][ jcoor ] * etaMass3/**/[ icoor ][ jcoor ][i][kcoor];    // \grad u^n:[/*2*/ * I\div d - (\grad d)^T] at each quadrature point
                    }
                auxMass3[ ig ][i][kcoor] = s;
                //if fullImplicit un=uk
            }

            ///////////////////////////////////////////////////////////////////////
        }

    }



    //
    // Numerical integration
    //
    Real g = 0.;
    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0; icoor < fe.nbLocalCoor(); icoor++ )
    {
        for ( kcoor = 0; kcoor < fe.nbLocalCoor(); kcoor++ )
        {

            // the block iccor of the elementary vector
            MatrixElemental::matrix_view mat = elmat.block ( icoor, kcoor );
            boost::shared_ptr<MatrixElemental::matrix_view> mat_convect;

            if (elmat_convect.get() )
            {
                mat_convect.reset (new MatrixElemental::matrix_view (elmat_convect->block ( icoor, kcoor ) ) );
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
                        for ( jcoor = 0; jcoor < fe.nbLocalCoor(); jcoor++ )
                        {
                            //source_mass1
                            s += convect_eta[ ig ][ jcoor ][j][kcoor] * guk[ ig ][ icoor ][jcoor] * fe.phi ( i, ig ) * fe.weightDet ( ig ) * rho;
                            //source_stress
                            s += B[ ig ][ icoor ][ jcoor ][j][kcoor] * fe.phiDer ( i, jcoor, ig ) * fe.weightDet ( ig );
                            //source_stress2
                            s -= fe.phiDer ( i, jcoor, ig ) * A2[ ig ][ icoor ][ jcoor ][j][kcoor] * fe.weightDet ( ig ) * mu;


                        }
                        //source_mass1
                        s += auxMass[ ig ][j][kcoor] * uk[ ig ][ icoor ] * fe.phi ( i, ig ) * fe.weightDet ( ig ) * rho;

                        //source_mass3
                        s += 0.5 * auxMass3[ ig ][j][kcoor] * uk[ ig ][ icoor ] * fe.phi ( i, ig ) * fe.weightDet ( ig ) * rho;

                        if (wImplicit)
                        {
                            //source_mass2{
                            s -= aux[ ig ][ icoor ][j][kcoor] * fe.phi ( i, ig ) * fe.weightDet ( ig ) * rho * alpha;
                            //convective term
                            g += aux[ ig ][ icoor ][j][kcoor] * fe.phi ( i, ig ) * fe.weightDet ( ig ) * rho;
                        }
                    }
                    mat ( i , j ) += s;
                    if (wImplicit && mat_convect.get() )
                    {
                        (*mat_convect) ( i , j ) += g;
                    }
                }
            }
        }
    }
}


//----------------------------------------------------------------------
// Compute the gradient in the Hdiv space, i.e. the opposite and transpose of the divergence matrix.
void grad_Hdiv ( Real coef, MatrixElemental& elmat, const CurrentFE& dualFE, const CurrentFE& primalFE, int iblock, int jblock )
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    Real sumdivphi (0.);
    const QuadratureRule& dualQuadRule ( dualFE.quadRule() );

    // Loop over all the degrees of freedom of the dual variable.
    for ( UInt i (0); i < dualFE.nbFEDof(); ++i )
    {
        // Loop over all the degrees of freedom of the primal variable.
        for ( UInt j (0); j < primalFE.nbFEDof(); ++j )
        {
            sumdivphi = 0.;
            // Loop over all the quadrature points.
            for ( UInt ig (0); ig < dualFE.nbQuadPt(); ++ig )
            {
                // There is no jacobian because of property of Piola transformation.
                sumdivphi -= primalFE.phi ( j, ig ) * dualFE.divPhiRef ( i, ig ) * dualQuadRule.weight ( ig );
            }
            mat ( i, j ) += coef * sumdivphi;
        }
    }
}

//----------------------------------------------------------------------
// Compute the divergence in the Hdiv space.
void div_Hdiv ( Real coef, MatrixElemental& elmat, const CurrentFE& dualFE, const CurrentFE& primalFE, int iblock, int jblock )
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    Real sumdivphi (0.);
    const QuadratureRule& dualQuadRule ( dualFE.quadRule() );

    // Loop over all the degrees of freedom of the dual variable.
    for ( UInt i (0); i < dualFE.nbFEDof(); ++i )
    {
        // Loop over all the degrees of freedom of the primal variable.
        for ( UInt j (0); j < primalFE.nbFEDof(); ++j )
        {
            sumdivphi = 0.;
            // Loop over all the quadrature points.
            for ( UInt ig (0); ig < dualFE.nbQuadPt(); ++ig )
            {
                // There is no jacobian because of property of Piola transformation.
                sumdivphi += primalFE.phi ( j, ig ) * dualFE.divPhiRef ( i, ig ) * dualQuadRule.weight ( ig );
            }
            mat ( j, i ) += coef * sumdivphi;
        }
    }
}

//----------------------------------------------------------------------
// Compute a Hdiv function dot product with the outwart unit normal times a hybrid function.
void TP_VdotN_Hdiv ( Real coef, MatrixElemental& elmat, const ReferenceFEHybrid& hybridFE,
                     const ReferenceFEHybrid& dualDotNFE, int iblock, int jblock )
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    UInt nbnode;
    Real sum (0.);

    // Loop over all the staticBdFE.
    for ( UInt nf (0); nf < hybridFE.numberBoundaryFE(); ++nf )
    {
        // Take the staticBdFE of the hybrid finite element.
        const CurrentFEManifold& boundaryElementHybridFE ( hybridFE[ nf ] );
        // Take the staticBdFE of the Hdiv function dot product outward unit normal.
        const CurrentFEManifold& boundaryElementDualDotNFE ( dualDotNFE[ nf ] );
        nbnode = boundaryElementHybridFE.nbFEDof();

        // Loop over all the the degrees of freedom of the dual dot normal variable.
        for ( UInt i (0); i < nbnode; ++i )
        {
            // Loop over all the degrees of freedom of the hybrid variable.
            for ( UInt j (0); j < nbnode; ++j )
            {
                sum = 0.;
                // Loop over all the quadrature point.
                for ( UInt ig (0); ig < boundaryElementHybridFE.nbQuadPt(); ++ig )
                {
                    // Using the Piola transform properties.
                    sum += boundaryElementHybridFE.phi ( j , ig ) *
                           boundaryElementDualDotNFE.phi ( i , ig ) *
                           boundaryElementHybridFE.wRootDetMetric ( ig );
                }
                // The matrix is block diagonal, so the size of the blocks is bdfe.nbNode.
                mat ( nf * nbnode + i, nf * nbnode + j ) += sum * coef;
            }
        }
    }
}

//----------------------------------------------------------------------
// Compute the mass matrix for hybrid variable.
void TP_TP_Hdiv ( Real coef, MatrixElemental& elmat, const ReferenceFEHybrid& hybridFE, int iblock, int jblock )
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    UInt nbnode;
    Real sum (0.);

    // Loop over all the staticBdFE.
    for ( UInt nf (0); nf < hybridFE.numberBoundaryFE(); ++nf )
    {
        // Take the staticBdFE of the hybrid finite element.
        const CurrentFEManifold& boundaryElementHybridFE ( hybridFE[ nf ] );
        nbnode = boundaryElementHybridFE.nbFEDof();

        // Loop over all the degrees of freedom of the first hybrid variable.
        for ( UInt i (0); i < nbnode; ++i )
        {
            // Loop over all the the degrees of freedom of the second hybrid variable.
            for ( UInt j (0); j < nbnode; ++j )
            {
                sum = 0.;
                // Loop over all the quadrature point.
                for ( UInt ig (0); ig < boundaryElementHybridFE.nbQuadPt() ; ++ig )
                    // Using the Piola transform properties.
                    sum += boundaryElementHybridFE.phi ( j , ig ) *
                           boundaryElementHybridFE.phi ( i , ig ) *
                           boundaryElementHybridFE.wRootDetMetric ( ig );

                // The matrix is block diagonal, so the size of the blocks is bdfe.nbNode.
                mat ( nf * nbnode + i, nf * nbnode + j ) += sum * coef;
            }
        }
    }
}


//----------------------------------------------------------------------
// Compute the mass matrix in Hdiv with a real scalar coefficient.
void mass_Hdiv ( Real coef, MatrixElemental& elmat, const CurrentFE& dualFE, int iblock, int jblock )
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    Real sum (0.);

    // Loop over all the degrees of freedom of the first dual variable.
    for ( UInt j (0); j < dualFE.nbFEDof(); ++j )
    {
        // Loop over all the degrees of freedom of the second dual variable.
        for ( UInt i (0); i < dualFE.nbFEDof() /* by symmetry j+1 */; ++i )
        {
            sum = 0.;
            // Loop over all the quadrature point.
            for ( UInt ig (0); ig < dualFE.nbQuadPt(); ++ig )
            {
                // Loop over all the space dimension, e.g. in 3D three times, in 2D two times.
                for ( UInt icoor (0); icoor < dualFE.nbLocalCoor() ; ++icoor )
                {
                    sum += dualFE.phi ( j, icoor, ig ) * dualFE.phi ( i, icoor, ig ) * dualFE.wDetJacobian ( ig );
                }
            }
            // Beware coef is the inverse of the permeability, not the permeability.
            mat ( i, j ) += sum * coef;
        }
    }
}

//----------------------------------------------------------------------
// Compute the mass matrix in Hdiv with a real matrix coefficient.
void mass_Hdiv ( Matrix const&  Invperm, MatrixElemental& elmat, const CurrentFE& dualFE, int iblock, int jblock )
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    Real partialSum (0.), sum (0.);

    // Loop over all the degrees of freedom of the first dual variable.
    for ( UInt j (0); j < dualFE.nbFEDof(); ++j )
    {
        // Loop over all the degrees of freedom of the second dual variable.
        for ( UInt i (0); i < dualFE.nbFEDof() /* by symmetry j+1 */; ++i )
        {
            sum = 0.;
            // Loop over all the quadrature point.
            for ( UInt ig (0); ig < dualFE.nbQuadPt(); ++ig )
            {
                // partialSum = phi[icoor]^T * K^-1 * phi[jcoor] for all icorr and jcoor.
                partialSum = 0.;
                /* Loop over all the space dimension of the first dual variable,
                    e.g. in 3D three times, in 2D two times.*/
                for ( UInt icoor (0); icoor < dualFE.nbLocalCoor(); ++icoor )
                {
                    /* Loop over all the space dimension of the first dual variable,
                          e.g. in 3D three times, in 2D two times.*/
                    for ( UInt jcoor (0); jcoor < dualFE.nbLocalCoor(); ++jcoor )
                    {
                        // Invperm is the inverse of the permeability.
                        partialSum += ( Invperm ( icoor, jcoor ) *
                                        dualFE.phi ( j, jcoor, ig ) *
                                        dualFE.phi ( i, icoor, ig ) );
                    }
                }
                sum += partialSum * dualFE.wDetJacobian ( ig );
            }
            mat ( i, j ) += sum ;
        }
    }
}

//----------------------------------------------------------------------
// Compute the mass matrix in Hdiv with a real function coefficient.
void mass_Hdiv ( Real ( *InvpermFun ) ( const Real&, const Real&, const Real& ),
                 MatrixElemental& elmat, const CurrentFE& dualFE, int iblock, int jblock )
{
    MatrixElemental::matrix_view mat = elmat.block ( iblock, jblock );
    Real sum (0.), x (0.), y (0.), z (0.);

    // Loop over all the degrees of freedom of the first dual variable.
    for ( UInt j (0); j < dualFE.nbFEDof(); ++j )
    {
        // Loop over all the degrees of freedom of the second dual variable.
        for ( UInt i (0); i < dualFE.nbFEDof() /* by symmetry j+1 */; ++i )
        {
            sum = 0.;
            // Loop over all the quadrature point.
            for ( UInt ig (0); ig < dualFE.nbQuadPt(); ++ig )
            {
                // Get the coordinate in the current element of the current quadrature point.
                x = dualFE.quadNode (ig, 0);
                y = dualFE.quadNode (ig, 1);
                z = dualFE.quadNode (ig, 2);

                // Loop over all the space dimension, e.g. in 3D three times, in 2D two times.
                for ( UInt icoor (0); icoor < dualFE.nbLocalCoor(); ++icoor )
                {
                    // Caution Invperm is the inverse of the permeability.
                    sum += InvpermFun ( x, y, z ) *
                           dualFE.phi ( j, icoor, ig ) *
                           dualFE.phi ( i, icoor, ig ) *
                           dualFE.wDetJacobian ( ig );
                }
            }
            mat ( i, j ) += sum;
        }
    }
}

//----------------------------------------------------------------------
// Compute the source vector in Hdiv with a real vector.
void source_Hdiv ( const Vector& source, VectorElemental& elvec, const CurrentFE& dualFE, int iblock )
{
    VectorElemental::vector_view vec = elvec.block ( iblock );
    Real sum (0.);

    // Loop over all the degrees of freedom of the dual variable.
    for ( UInt i (0); i < dualFE.nbFEDof(); ++i )
    {
        // Loop over all the quadrature point.
        for ( UInt ig (0); ig < dualFE.nbQuadPt(); ++ig )
        {
            sum = 0.;
            /* Loop over all the space dimension of the dual variable,
               e.g. in 3D three times, in 2D two times.*/
            for ( UInt icoor (0); icoor < dualFE.nbLocalCoor(); ++icoor )
            {
                // sum[icoor] = source[icoor] * phi[icoor] for all icorr.
                sum += source ( icoor ) * dualFE.phi ( i, icoor, ig );

            }

            vec ( i ) += sum * dualFE.wDetJacobian ( ig );
        }
    }
}

} // namespace LifeV

#endif
