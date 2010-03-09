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

/*!
    @file elemOperHdiv.cpp
    @brief functions for Hdiv that were in elemOper.hpp

    @author Simone Deparis <simone.deparis@epfl.ch>
    @date 09 Mar 2010

 */

#include <life/lifefem/elemOperHdiv.hpp>

namespace LifeV {

//----------------------------------------------------------------------
/*! \function grad_Hdiv : compute
  - coef * \int_{current element} q_j * div w_i
  where w_j is a vectorial H(div) basis function,
  and q_j is a L2 basis function.
  \param coef  : constant coefficient.
  \param elmat : (mixed) element matrix.
  \param fe_u  : current vectorial element (in H(div))
  \param fe_p  : current scalar element (in L2)
  \param iblock, \param jblock : subarray indexes where to store the integral just computed.
*/
void grad_Hdiv( Real coef, ElemMat& elmat, const CurrentHdivFE& fe_u,
                const CurrentFE& fe_p, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int ig, i, j;
    Real sumdivphi;
    for ( i = 0;i < fe_u.nbNode;i++ )
    {
        for ( j = 0;j < fe_p.nbNode;j++ )
        {
            sumdivphi = 0;
            for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )
            {
                sumdivphi -= fe_p.phi( j, ig ) * fe_u.divPhi( i, ig ) * fe_u.qr.weight( ig );
                //! there is no jacobian because of property of Piola transformation
            }
            mat( i, j ) += coef * sumdivphi;
        }
    }
}

//----------------------------------------------------------------------
/*! \function div_Hdiv : compute
  coef * \int_{current element} q_j * div w_i
  where w_j is a vectorial H(div) basis function,
  and q_j is a scalar L2 basis function.
  \param coef  : constant coefficient.
  \param elmat : (mixed) element matrix.
  \param fe_u  : current vectorial element (in H(div))
  \param fe_p  : current scalar element (in L2)
  \param iblock, \param jblock : subarray indexes where to store the integral just computed.
*/
void div_Hdiv( Real coef, ElemMat& elmat, const CurrentHdivFE& fe_u,
               const CurrentFE& fe_p, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int ig, i, j;
    Real sumdivphi;
    for ( i = 0;i < fe_u.nbNode;i++ )
    {
        for ( j = 0;j < fe_p.nbNode;j++ )
        {
            sumdivphi = 0;
            for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )
            {
                sumdivphi += fe_p.phi( j, ig ) * fe_u.divPhi( i, ig ) * fe_u.qr.weight( ig );
                //! there is no jacobian because of property of Piola transformation
            }
            mat( j, i ) += coef * sumdivphi;
        }
    }
}

//----------------------------------------------------------------------
/*! \function TP_VdotN_Hdiv : compute

coef * \int_{BOUNDARY of current element} lambda_j * { w_i \cdot n }

where w_j is a vectorial H(div) basis function,
and lambda_j are the Lagrange multiplier basis functions
that enforce continuity of the normal component of
the vectorial functions across two neighbouring elements.
Interprated as trace of pressure... (TP)
See Hybridization for Mixed Hybrid Finite Element Method.

Thanks to the Piola transform, the computation is performed
on the boundary of the REFERENCE Element. But in general, the
boundary of a 3D Reference element is not a 2D Reference element.
Example:
REFERENCE TETRA -> 3 REFERENCE TRIA + 1 EQUILATERAL TRIANGLE...
REFERENCE PRISM -> 2 TRIA + 3 QUAD...?
REFERENCE HEXA  -> 6 REFERENCE QUAD.

n : is the normal unit vector oriented outward of the current element.

\param coef  : constant coefficient.
\param elmat : (mixed) element matrix.
\param tpfe  : reference lagrange multiplier element (for hybrid MFE)
\param iblock, \param jblock : subarray indexes where to store the integral just computed.
*/
void TP_VdotN_Hdiv( Real coef, ElemMat& elmat, const RefHybridFE& tpfe,
                    const RefHybridFE& vdotnfe, int iblock, int jblock )
{
    //! previous way of construction (worked only for RTO hexa)
    // TP_TP_Hdiv(coef, elmat, tpfe, iblock, jblock);

    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Int ig, i, j, nbnode;
    UInt nf;
    Real tpvn;
    // Real a[4] = {2./sqrt(3.),2,2,2};

    for ( nf = 0 ; nf < tpfe.nBdFE() ; nf ++ )
    {
        //! use the static boundary element of the reference element.
        const StaticBdFE & bdtpfe = tpfe[ nf ];
        const StaticBdFE & bdvdotnfe = vdotnfe[ nf ];
        nbnode = bdtpfe.nbNode;

        for ( i = 0 ; i < nbnode ; i ++ )
        {
            for ( j = 0 ; j < nbnode ; j ++ )
            {
                tpvn = 0.;
                for ( ig = 0 ; ig < bdtpfe.nbQuadPt ; ig ++ )
                    tpvn += /*a[nf]*/ bdtpfe.phi( j , ig ) * bdvdotnfe.phi( i , ig ) * bdtpfe.weightMeas( ig );
                //! using the Piola transform properties.

                //! Matrix : block diagonal. size of the blocks = bdfe.nbNode.
                mat( nf * nbnode + i, nf * nbnode + j ) += tpvn * coef;
            }
        }
    }
}

//----------------------------------------------------------------------
/*! \function TP_TP_Hdiv : compute

coef * \int_{BOUNDARY of current element} lambda_j * lambda_i

where lambda_j are the Lagrange multiplier basis functions
that enforce continuity of the normal component of
the vectorial functions across two neighbouring elements.
Interprated as trace of pressure... (TP)
See Hybridization for Mixed Hybrid Finite Element Method.

Thanks to the Piola transform, the computation is performed
on the boundary of the REFERENCE Element. But in general, the
boundary of a 3D Reference element is not a 2D Reference element.
Example:
REFERENCE TETRA -> 3 REFERENCE TRIA + 1 EQUILATERAL TRIANGLE...
REFERENCE PRISM -> 2 TRIA + 3 QUAD...?
REFERENCE HEXA  -> 6 REFERENCE QUAD.

n : is the normal unit vector oriented outward of the current element.

\param coef  : constant coefficient.
\param elmat : (mixed) element matrix.
\param tpfe  : reference lagrange multiplier element (for hybrid MFE)
\param iblock, \param jblock : subarray indexes where to store the integral just computed.
*/
void TP_TP_Hdiv( Real coef, ElemMat& elmat, const RefHybridFE& tpfe, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Int ig, i, j, nbnode;
    UInt nf;
    Real tpvn;

    for ( nf = 0 ; nf < tpfe.nBdFE() ; nf ++ )
    {
        //! use the static boundary element of the reference element.
        const StaticBdFE & bdfe = tpfe[ nf ];
        nbnode = bdfe.nbNode;

        for ( i = 0 ; i < nbnode ; i ++ )
        {
            for ( j = 0 ; j < nbnode ; j ++ )
            {
                tpvn = 0.;
                for ( ig = 0 ; ig < bdfe.nbQuadPt ; ig ++ )
                    tpvn += bdfe.phi( j , ig ) * bdfe.phi( i , ig ) * bdfe.weightMeas( ig );
                //! using the Piola transform properties.

                //! Matrix : block diagonal. size of the blocks = bdfe.nbNode.
                mat( nf * nbnode + i, nf * nbnode + j ) += tpvn * coef;
            }
        }
    }
}

//----------------------------------------------------------------------
/*! \function mass_Hdiv : compute

coef *  \int_{current element} w_j * w_i
where w_j is a vectorial H(div) basis function

Here the permeability matrix is a CONSTANT SCALAR tensor (i.e. = coef * Id).

BEWARE  :   it is "coef" that is used (and NOT ITS INVERSE!!).

\param coef  : constant coefficient.
\param elmat : (mixed) element matrix.
\param fe_u  : current vectorial element (in H(div))
\param iblock, \param jblock : subarray indexes where to store the integral just computed.
*/
void mass_Hdiv( Real coef, ElemMat& elmat, const CurrentHdivFE& fe, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int ig, i, j, icoor;
    Real x;
    for ( j = 0 ; j < fe.nbNode ; j ++ )
    {
        for ( i = 0 ; i < fe.nbNode /* by symmetry j+1 */ ; i ++ )
        {
            x = 0.;
            for ( ig = 0 ; ig < fe.nbQuadPt ; ig ++ )
            {
                for ( icoor = 0 ; icoor < fe.nbCoor ; icoor ++ )
                {
                    x += fe.phi( j , icoor , ig ) * fe.phi( i , icoor , ig ) * fe.weightDet( ig );
                }
            }
            //! coef is the inverse of the permeability
            mat( i, j ) += x * coef;
        }
    }
}

//----------------------------------------------------------------------
/*! \function mass_Hdiv : compute
  \int_{current element} ((Invperm* w_j) * w_i
  where w_j is a vectorial H(div) basis function

  Here the permeability matrix "Invperm" is a CONSTANT symmetric positive definite
  matrix (NON DIAGONAL a priori). The matrix is constant over the
  whole current element and is already inverted by LU or Choleski Lapack.

  \param Invperm : constant coefficient TENSOR. (CONSTANT over the current element).
  \param elmat   : (mixed) element matrix.
  \param fe      : current vectorial element (in H(div))
  \param iblock, \param jblock : subarray indexes where to store the integral just computed.
*/
void mass_Hdiv( Matrix const&  Invperm, ElemMat& elmat, const CurrentHdivFE& fe,
                int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    for (Int j = 0 ; j < fe.nbNode ; j ++ )
    {
        for (Int i = 0 ; i < fe.nbNode /* by symmetry j+1 */ ; i ++ )
        {
            Real x = 0.;
            for (Int ig = 0 ; ig < fe.nbQuadPt ; ig ++ )
            {
                double __t = 0;
                //x = phi[icoor]^T * K^-1 * phi[jcoor]
                for (Int icoor = 0 ; icoor < fe.nbCoor ; icoor ++ )
                {
                    for (Int jcoor = 0 ; jcoor < fe.nbCoor ; jcoor ++ )
                    {
                        //! Invperm is the inverse of the permeability
                        __t += (  Invperm( icoor, jcoor ) *
                                  fe.phi( j, jcoor, ig ) *
                                  fe.phi( i, icoor, ig ) );
                    }
                }
                x += __t*fe.weightDet( ig );
            }
            mat( i, j ) += x ;
        }
    }
}


//----------------------------------------------------------------------
/*! \function mass_Hdiv : compute
  \int_{current element} ((Invperm* w_j) * w_i
  where w_j is a vectorial H(div) basis function

  Here the permeability is a NON-CONSTANT scalar function
  whose inverse is "Invperm".

  We note again that it is the inverse of the permeability that
  is provided directly (Invperm = K^{-1}).

  \param Invperm : scalar function inverse of the permeability.
  \param elmat   : (mixed) element matrix.
  \param fe      : current vectorial element (in H(div))
  \param iblock, \param jblock : subarray indexes where to store the integral just computed.
*/
void mass_Hdiv( Real ( *Invperm ) ( const Real&, const Real&, const Real& ),
                ElemMat& elmat, const CurrentHdivFE& fe, int iblock, int jblock )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int ig, i, j, icoor;
    Real intg, x, y, z;
    for ( j = 0 ; j < fe.nbNode ; j ++ )
    {
        for ( i = 0 ; i < fe.nbNode /* by symmetry j+1 */ ; i ++ )
        {
            intg = 0.;
            for ( ig = 0 ; ig < fe.nbQuadPt ; ig ++ )
            {
                fe.coorQuadPt( x, y, z, ig );
                for ( icoor = 0 ; icoor < fe.nbCoor ; icoor ++ )
                {
                    //! caution inverse of the permeability
                    intg += 1. * Invperm( x, y, z ) *
                        fe.phi( j , icoor , ig ) *
                        fe.phi( i , icoor , ig ) *
                        fe.weightDet( ig );
                }
            }
            mat( i, j ) += intg;
        }
    }
}

/*
  \int_{current element} w_j * w_i where w_j is a H(div) basis function.

  Here the permeability matrix is a NON-CONSTANT symmetric positive definite
  matrix (NON DIAGONAL a priori).
*/

/*

TO BE DONE LATER... IN ELEMOPER.H

void mass_Hdiv(KNM<Real> &Kperm, ElemMat& elmat, const CurrentHdivFE& fe,
int iblock=0, int jblock=0)
{
ElemMat::matrix_view mat = elmat.block(iblock,jblock);
// int ig,i,j,icoor;
Real s;
KN<Real> p(fe.nbCoor);
KN<Real> b(fe.nbCoor);
KN<Real> x(fe.nbCoor);
KN<Real> y(fe.nbCoor);

for(int j=0;j<fe.nbNode;j++){
for(int i=0;i<fe.nbNode;i++){
s =0.;
for(int ig=0;ig<fe.nbQuadPt;ig++){
choldc(Kperm, p);
for(int lcoor=0;lcoor<fe.nbCoor;lcoor++){
b(lcoor) = fe.phi(j,lcoor,ig);
}
cholsl(Kperm,p, b, y);
for(int icoor=0;icoor<fe.nbCoor;icoor++){
s += y(icoor)*fe.phi(i,icoor,ig)*fe.weightDet(ig);
}
}
mat(i,j) += s;
}
}
}
*/


void mass_Mixed_Hdiv( Real coef, ElemMat& elmat, const CurrentFE& fe,
                      const CurrentHdivFE& hdivfe, int iblock, int jblock )
{
    // to be improved: symmetry not used
    Real x;
    for ( int icoor = 0;icoor < fe.nbCoor; icoor ++ )
    {
        ElemMat::matrix_view mat = elmat.block( iblock + icoor, jblock );
        for ( int j = 0 ; j < hdivfe.nbNode ; j ++ )
        {
            for ( int i = 0 ; i < fe.nbNode ; i ++ )
            {
                x = 0.;
                for ( int ig = 0 ; ig < fe.nbQuadPt ; ig ++ )
                {
                    x += hdivfe.phi( j , icoor , ig ) * fe.phi( i , ig ) * fe.weightDet( ig );
                }
                mat( i, j ) += coef * x;
            }
        }
    }
}

// ===================================================
// Constructors & Destructor
// ===================================================
elemOperHdiv::elemOperHdiv() :
    M_variableOne (),
    M_variableTwo ()
{

}

elemOperHdiv::elemOperHdiv( first_Type&  variableOne,
                              second_Type& variableTwo ) :
    M_variableOne ( variableOne ),
    M_variableTwo ( variableTwo )

} // Namespace LifeV
