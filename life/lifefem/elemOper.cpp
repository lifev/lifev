/*-*- mode: c++ -*-
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!@file elemOper.cpp
   Element matrix operations
*/

#include <life/lifefem/elemOper.hpp>

namespace LifeV
{
//
//----------------------------------------------------------------------
//                      Element matrix operator
//----------------------------------------------------------------------
//Real coef(t,x,y,z,u)
/*
  Mass matrix: \int coef(t,x,y,z,u) v_i v_j
*/
void mass( Real (*coef)(Real,Real,Real,Real,Real),
           ElemMat& elmat, const CurrentFE& fe,
           const Dof& dof,
           const ScalUnknown<Vector>& U,Real t)
{
    ASSERT_PRE( fe.hasJac(), "Mass matrix needs at least the jacobian" );
    int iblock=0,jblock=0;
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int i, ig;
    int iloc, jloc;
    Real s, coef_s;
    ID eleId=fe.currentId();
    int iu;
    double uPt;
    Real x,y,z;

    std::vector<Real> locU(fe.nbNode);
    for (i=0;i<fe.nbNode;i++)
    {
        locU[i]=U[dof.localToGlobal(eleId,i+1)-1];    //(one component)
    }

    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            uPt=0.0;
            for(iu=0;iu<fe.nbNode;iu++){
                uPt+=locU[iu]*fe.phi(iu,ig);
            }
            fe.coorQuadPt(x,y,z,ig);
            s += fe.phi( iloc, ig ) * fe.phi( iloc, ig ) * fe.weightDet( ig )*
                coef(t,x,y,z,uPt);
        }
        mat( iloc, iloc ) += s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            uPt=0.0;
            for(iu=0;iu<fe.nbNode;iu++){
                uPt+=locU[iu]*fe.phi(iu,ig);
            }
            fe.coorQuadPt(x,y,z,ig);
            s += fe.phi( iloc, ig ) * fe.phi( jloc, ig ) * fe.weightDet( ig )*
                coef(t,x,y,z,uPt);
        }
        coef_s = s;
        mat( iloc, jloc ) += coef_s;
        mat( jloc, iloc ) += coef_s;
    }
}


/*
  Stiffness matrix: \int coef(t,x,y,z,u) grad v_i . grad v_j
*/
void stiff( Real (*coef)(Real,Real,Real,Real,Real),
            ElemMat& elmat, const CurrentFE& fe,
            const Dof& dof,
            const ScalUnknown<Vector>& U,Real t)
{
    int iblock=0,jblock=0;
    ASSERT_PRE( fe.hasFirstDeriv(),
                "Stiffness matrix needs at least the first derivatives" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int iloc, jloc;
    int i, icoor, ig;
    double s, coef_s;
    ID eleId=fe.currentId();
    int iu;
    double uPt;
    Real x,y,z;

    std::vector<Real> locU(fe.nbNode);
    for (i=0;i<fe.nbNode;i++)
    {
        locU[i]=U[dof.localToGlobal(eleId,i+1)-1];    //(one component)
    }
    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            uPt=0.0;
            for(iu=0;iu<fe.nbNode;iu++){
                uPt+=locU[iu]*fe.phi(iu,ig);
            }
            fe.coorQuadPt(x,y,z,ig);
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig )
                    * fe.weightDet( ig )*coef(t,x,y,z,uPt);
        }
        mat( iloc, iloc ) += s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            uPt=0.0;
            for(iu=0;iu<fe.nbNode;iu++){
                uPt+=locU[iu]*fe.phi(iu,ig);
            }
            fe.coorQuadPt(x,y,z,ig);
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                    fe.weightDet( ig )*coef(t,x,y,z,uPt);
        }
        coef_s = s;
        mat( iloc, jloc ) += coef_s;
        mat( jloc, iloc ) += coef_s;
    }
}


/*
  compute \int fct(t,x,y,z,u) \phi_i
*/
void source( Real (*fct)(Real,Real,Real,Real,Real),
             ElemVec& elvec, const CurrentFE& fe,
             const Dof& dof,
             const ScalUnknown<Vector>& U,Real t)
{
    int iblock=0;
    ASSERT_PRE( fe.hasQuadPtCoor(),
                "Source with space dependent function need updated quadrature "
                "point coordinates. Call for example updateFirstDerivQuadPt() "
                "instead of updateFirstDeriv()." );
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;
    ID eleId=fe.currentId();
    int iu;
    double uPt;

    std::vector<Real> locU(fe.nbNode);
    for (i=0;i<fe.nbNode;i++)
    {
        locU[i]=U[dof.localToGlobal(eleId,i+1)-1];    //(one component)
    }
    for ( i = 0;i < fe.nbNode;i++ )
    {
        s = 0.0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            uPt=0.0;
            for(iu=0;iu<fe.nbNode;iu++){
                uPt+=locU[iu]*fe.phi(iu,ig);
            }
            s += fe.phi( i, ig ) *
                fct(t, fe.quadPt( ig, 0 ),
                    fe.quadPt( ig, 1 ),
                    fe.quadPt( ig, 2 ), uPt) *
                fe.weightDet( ig );
        }
        vec( i ) += s;
    }
}











//
// coeff*Mass
//
void mass( Real coef, ElemMat& elmat, const CurrentFE& fe,
           int iblock, int jblock )
/*
  Mass matrix: \int v_i v_j
*/
{
    ASSERT_PRE( fe.hasJac(), "Mass matrix needs at least the jacobian" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int i, ig;
    int iloc, jloc;
    Real s, coef_s;
    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s += fe.phi( iloc, ig ) * fe.phi( iloc, ig ) * fe.weightDet( ig );
        }
        mat( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
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
           int iblock, int jblock, int nb )
/*
  Mass matrix: \int v_i v_j (nb blocks on the diagonal, nb>1)
*/
{
    ASSERT_PRE( fe.hasJac(), "Mass matrix needs at least the jacobian" );
    ASSERT_PRE( nb > 1, "if nb = 1, use the other mass function" );
    Tab2d mat_tmp( fe.nbNode, fe.nbNode );
    int i, ig;
    int iloc, jloc;
    Real s, coef_s;
    mat_tmp = ZeroMatrix( fe.nbNode, fe.nbNode );
    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s += fe.phi( iloc, ig ) * fe.phi( iloc, ig ) * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
            s += fe.phi( iloc, ig ) * fe.phi( jloc, ig ) * fe.weightDet( ig );
        coef_s = coef * s;
        mat_tmp( iloc, jloc ) += coef_s;
        mat_tmp( jloc, iloc ) += coef_s;
    }
    // copy on the components
    for ( int icomp = 0;icomp < nb;icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0;i < fe.nbDiag;i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
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
           int iblock, int jblock, int nb )
/*
  Mass matrix: \int v_i v_j (nb blocks on the diagonal, nb>1)
*/
{
    ASSERT_PRE( fe.hasJac(), "Mass matrix needs at least the jacobian" );
    ASSERT_PRE( nb > 1, "if nb = 1, use the other mass function" );
    Tab2d mat_tmp( fe.nbNode, fe.nbNode );
    int i, ig;
    int iloc, jloc;
    Real s;//, coef_s;
    mat_tmp = ZeroMatrix( fe.nbNode, fe.nbNode );
    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s += coef[ig] * fe.phi( iloc, ig ) * fe.phi( iloc, ig ) * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) = s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
            s += coef[ig]*fe.phi( iloc, ig ) * fe.phi( jloc, ig ) * fe.weightDet( ig );
        mat_tmp( iloc, jloc ) += s;
        mat_tmp( jloc, iloc ) += s;
    }
    // copy on the components
    for ( int icomp = 0;icomp < nb;icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0;i < fe.nbDiag;i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
        {
            iloc = fe.patternFirst( i );
            jloc = fe.patternSecond( i );
            mat_icomp( iloc, jloc ) += mat_tmp( iloc, jloc );
            mat_icomp( jloc, iloc ) += mat_tmp( jloc, iloc );
        }
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

    ASSERT_PRE( fe1.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );
    ASSERT_PRE( fe2.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );

    ElemMat::matrix_view mat = elmat.block( iblock, jblock );



    Real sum, sum1, sum2;
    int i, j, ig, icoor, jcoor;
    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];
    Real phid1[ fe1.nbNode ][ fe1.nbCoor ][ bdfe.nbQuadPt ];
    Real phid2[ fe2.nbNode ][ fe2.nbCoor ][ bdfe.nbQuadPt ];
    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( int ig = 0; ig < bdfe.nbQuadPt; ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbNode; ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                drp1[ icoor ] = fe1.refFE.dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE.dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
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
    for ( i = 0; i < fe1.nbNode; ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbNode; ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt ; ++ig )
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


    ASSERT_PRE( fe1.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );
    ASSERT_PRE( fe2.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );

    ElemMat::matrix_type mat_tmp( fe1.nbNode, fe2.nbNode );


    Real sum, sum1, sum2;
    int i, j, ig, icoor, jcoor;
    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];
    Real phid1[ fe1.nbNode ][ fe1.nbCoor ][ bdfe.nbQuadPt ];
    Real phid2[ fe2.nbNode ][ fe2.nbCoor ][ bdfe.nbQuadPt ];
    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( int ig = 0; ig < bdfe.nbQuadPt; ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbNode; ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                drp1[ icoor ] = fe1.refFE.dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE.dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
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
    for ( i = 0; i < fe1.nbNode; ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbNode; ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt ; ++ig )
                    sum += phid1[ i ][ icoor ][ ig ] * phid2[ j ][ icoor ][ ig ] * bdfe.weightMeas( ig );
            mat_tmp( i, j ) = coef * sum;
        }
    }


    // copy on the components
    for ( int icomp = 0;icomp < nb;icomp++ )
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

    ASSERT_PRE( fe1.hasFirstDeriv(),
                "ipstab_bgrad needs at least the first derivatives" );
    ASSERT_PRE( fe2.hasFirstDeriv(),
                "ipstab_bgrad needs at least the first derivatives" );

    ElemMat::matrix_type mat_tmp( fe1.nbNode, fe2.nbNode );

    Real sum, sum1, sum2;
    int i, j, icoor, jcoor, ig;

    //
    // convection velocity \beta on the boundary quadrature points
    //
    Real b[ fe1.nbCoor ][ bdfe.nbQuadPt ];

    for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
    {
        for ( ig = 0;ig < bdfe.nbQuadPt;ig++ )
        {
            sum = 0;
            for ( i = 0; i < bdfe.nbNode; ++i )
            {
                sum += bdfe.phi( i, ig ) * beta.vec() [ icoor * bdfe.nbNode + i ];
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
    Real phid1[ fe1.nbNode ][ fe1.nbCoor ][ bdfe.nbQuadPt ];
    Real phid2[ fe2.nbNode ][ fe2.nbCoor ][ bdfe.nbQuadPt ];
    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( int ig = 0; ig < bdfe.nbQuadPt; ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbNode; ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                drp1[ icoor ] = fe1.refFE.dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE.dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
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
    for ( i = 0; i < fe1.nbNode; ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbNode; ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
                for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
                    for ( ig = 0;ig < bdfe.nbQuadPt;ig++ )
                        sum += phid1[ i ][ icoor ][ ig ]*phid2[ j ][ jcoor ][ ig ]
                            *b[ icoor ][ ig ]*b[ jcoor ][ ig ]
                            *bdfe.weightMeas( ig );
            mat_tmp( i, j ) = coef * sum;
        }
    }

    // copy on the components
    for ( int icomp = 0;icomp < nb; icomp++ )
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

    ASSERT_PRE( fe1.hasFirstDeriv(),
                "ipstab_div needs at least the first derivatives" );
    ASSERT_PRE( fe2.hasFirstDeriv(),
                "ipstab_div needs at least the first derivatives" );


    Real sum, sum1, sum2;
    int i, j, ig, icoor, jcoor;
    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];
    Real phid1[ fe1.nbNode ][ fe1.nbCoor ][ bdfe.nbQuadPt ];
    Real phid2[ fe2.nbNode ][ fe2.nbCoor ][ bdfe.nbQuadPt ];
    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( int ig = 0; ig < bdfe.nbQuadPt; ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbNode; ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                drp1[ icoor ] = fe1.refFE.dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE.dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
                {
                    sum1 += fe1.tInvJac( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac( icoor, jcoor, 0 ) * drp2[ jcoor ];
                }
                phid1[ i ][ icoor ][ ig ] = sum1;
                phid2[ i ][ icoor ][ ig ] = sum2;
            }
        }
    }

    for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
    {
        for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
        {
            ElemMat::matrix_view mat_icomp = elmat.block( iblock + icoor, jblock + jcoor );
            // Loop on rows
            for ( i = 0; i < fe1.nbNode; ++i )
            {
                // Loop on columns
                for ( j = 0; j < fe2.nbNode; ++j )
                {
                    sum = 0.0;
                    for ( ig = 0;ig < bdfe.nbQuadPt;ig++ )
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

    ASSERT_PRE( fe1.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );
    ASSERT_PRE( fe2.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );

    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    Real sum, sum1, sum2;
    int i, j, ig, icoor, jcoor;
    Real phid1[ fe1.nbNode ][ fe1.nbCoor ][ bdfe.nbQuadPt ];
    Real phid2[ fe2.nbNode ][ fe2.nbCoor ][ bdfe.nbQuadPt ];

    std::vector<Real> x(3), rx1(3), drp1(3), rx2(3), drp2(3);
    std::vector<Real> b1(3), b2(3);

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2
    const KNM<Real>& normal = bdfe.normal;

    //
    // convection velocity term |\beta . n|^2 / |\beta|
    // on the boundary quadrature points
    //
    Real ba2[ bdfe.nbQuadPt ];

    for ( ig = 0; ig < bdfe.nbQuadPt; ig++ )
    {
        sum1 = 0;
        sum2 = 0;
        for ( icoor = 0; icoor < fe1.nbCoor; ++icoor ) {
            for ( i = 0; i < bdfe.nbNode; ++i ) {
                Real betaLoc = bdfe.phi( i, ig ) *
                    beta.vec() [ icoor * bdfe.nbNode + i ];
                sum1 += betaLoc * normal(icoor, ig);
                sum2 += betaLoc * betaLoc;
            }
        }
        ba2[ ig ] = sum2 == 0 ? 0 : sum1 * sum1 / pow( sum2, 0.5 );
    }

    for ( int ig = 0; ig < bdfe.nbQuadPt; ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbNode; ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                drp1[ icoor ] = fe1.refFE.dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE.dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
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
    for ( i = 0; i < fe1.nbNode; ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbNode; ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt ; ++ig )
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

    ASSERT_PRE( fe1.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );
    ASSERT_PRE( fe2.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );

    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    Real sum, sum1, sum2;
    int i, j, ig, icoor, jcoor;
    Real phid1[ fe1.nbNode ][ fe1.nbCoor ][ bdfe.nbQuadPt ];
    Real phid2[ fe2.nbNode ][ fe2.nbCoor ][ bdfe.nbQuadPt ];

    std::vector<Real> x(3), rx1(3), drp1(3), rx2(3), drp2(3);
    std::vector<Real> b1(3), b2(3);

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2
    const KNM<Real>& normal = bdfe.normal;

    //
    // convection velocity term |\beta . n|
    // on the boundary quadrature points
    //
    Real bn[ bdfe.nbQuadPt ];

    for ( ig = 0; ig < bdfe.nbQuadPt; ig++ )
    {
        sum1 = 0;
        sum2 = 0;

        for ( icoor = 0; icoor < fe3.nbCoor; ++icoor ) {
            for ( i = 0; i < fe3.nbNode; ++i ) {
	        Real betaLoc = fe3.phi( i, ig ) *
		  beta.vec() [ icoor*fe3.nbNode + i ];
                sum1 += betaLoc * normal(icoor, ig);
            }
        }
        bn[ ig ] = std::abs(sum1);
    }

    for ( int ig = 0; ig < bdfe.nbQuadPt; ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbNode; ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                drp1[ icoor ] = fe1.refFE.dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE.dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
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
    for ( i = 0; i < fe1.nbNode; ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbNode; ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt ; ++ig )
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
    ASSERT_PRE( fe.hasFirstDeriv(),
                "Stiffness matrix needs at least the first derivatives" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int iloc, jloc;
    int i, icoor, ig;
    double s, coef_s;
    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig )
                    * fe.weightDet( ig );
        }
        mat( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
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
    ASSERT_PRE( fe.hasFirstDeriv() && fe.hasQuadPtCoor(),
                "Stiffness matrix with a diffusion function needs the first derivatives and the coordinates of the quadrature points.  Call for example updateFirstDerivQuadPt() instead of updateFirstDeriv()" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int iloc, jloc;
    int i, icoor, ig;
    double s, coef_s, coef_f;
    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            coef_f = fct( fe.quadPt( ig, 0 ), fe.quadPt( ig, 1 ), fe.quadPt( ig, 2 ) );
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                s += coef_f * fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig )
                    * fe.weightDet( ig );
        }
        mat( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            coef_f = fct( fe.quadPt( ig, 0 ), fe.quadPt( ig, 1 ), fe.quadPt( ig, 2 ) );
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
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


    ASSERT_PRE( fe.hasFirstDeriv(),
                "Stiffness (vect) matrix needs at least the first derivatives" );
    ASSERT_PRE( nb > 1, "if nb = 1, use the other stiff function" );

    Tab2d mat_tmp( fe.nbNode, fe.nbNode );
    mat_tmp = ZeroMatrix( fe.nbNode, fe.nbNode );

    int iloc, jloc;
    int i, icoor, ig;
    double s, coef_s;
    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig )
                    * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                    fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp( iloc, jloc ) += coef_s;
        mat_tmp( jloc, iloc ) += coef_s;
    }
    // copy on the components
    for ( int icomp = 0;icomp < nb;icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0;i < fe.nbDiag;i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
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

    ASSERT_PRE( fe.hasFirstDeriv(),
                "Stiffness (vect) matrix needs at least the first derivatives" );
    ASSERT_PRE( nb > 1, "if nb = 1, use the other stiff function" );

    Tab2d mat_tmp( fe.nbNode, fe.nbNode );
    mat_tmp = ZeroMatrix( fe.nbNode, fe.nbNode );

    int iloc, jloc;
    int i, icoor, ig;
    double s;//, coef_s;
    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                s += coef[ig] *fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig )
                    * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) += s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                s += coef[ig] * fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                    fe.weightDet( ig );
        }
        mat_tmp( iloc, jloc ) += s;
        mat_tmp( jloc, iloc ) += s;
    }
    // copy on the components
    for ( int icomp = 0;icomp < nb;icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0;i < fe.nbDiag;i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
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
                 int iblock, int jblock, int nb )


/*
  Stiffness matrix: coef*\int curl v_i . curl v_j (nb blocks on the diagonal, nb>1)
*/
{


    ASSERT_PRE( fe.hasFirstDeriv(),
                "Stiffness (vect) matrix needs at least the first derivatives" );
    ASSERT_PRE( nb > 1, "if nb = 1, use the other stiff function" );
    Tab2d mat_tmp( fe.nbNode, fe.nbNode );
    mat_tmp = ZeroMatrix( fe.nbNode, fe.nbNode );
    Tab2d mat_tmp11( fe.nbNode, fe.nbNode );
    mat_tmp11 = ZeroMatrix( fe.nbNode, fe.nbNode );
    Tab2d mat_tmp12( fe.nbNode, fe.nbNode );
    mat_tmp12 = ZeroMatrix( fe.nbNode, fe.nbNode );
    Tab2d mat_tmp13( fe.nbNode, fe.nbNode );
    mat_tmp13 = ZeroMatrix( fe.nbNode, fe.nbNode );
    Tab2d mat_tmp21( fe.nbNode, fe.nbNode );
    mat_tmp21 = ZeroMatrix( fe.nbNode, fe.nbNode );
    Tab2d mat_tmp22( fe.nbNode, fe.nbNode );
    mat_tmp22 = ZeroMatrix( fe.nbNode, fe.nbNode );
    Tab2d mat_tmp23( fe.nbNode, fe.nbNode );
    mat_tmp23 = ZeroMatrix( fe.nbNode, fe.nbNode );
    Tab2d mat_tmp31( fe.nbNode, fe.nbNode );
    mat_tmp31 = ZeroMatrix( fe.nbNode, fe.nbNode );
    Tab2d mat_tmp32( fe.nbNode, fe.nbNode );
    mat_tmp32 = ZeroMatrix( fe.nbNode, fe.nbNode );
    Tab2d mat_tmp33( fe.nbNode, fe.nbNode );
    mat_tmp33 = ZeroMatrix( fe.nbNode, fe.nbNode );



    int iloc, jloc;
    int i, ig;
    double s, coef_s;

    // diagonal 11
    //
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt; ig++ )
        {
            s = fe.phiDer( iloc, 1, ig ) * fe.phiDer( iloc, 1, ig ) * fe.weightDet( ig )
                + fe.phiDer( iloc, 2, ig ) * fe.phiDer( iloc, 2, ig ) * fe.weightDet( ig ) ;
        }
        mat_tmp11( iloc, iloc ) += coef * s;
    }
    // extra diagonal 11
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
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
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt; ig++ )
        {
            s = - fe.phiDer( iloc, 1, ig ) * fe.phiDer( iloc, 0, ig ) * fe.weightDet( ig );
        }
        mat_tmp12( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 12
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s = - fe.phiDer( iloc, 1, ig ) * fe.phiDer( jloc, 0, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp12( iloc, jloc ) -= coef_s;
        mat_tmp12( jloc, iloc ) -= coef_s;
    }

    // diagonal 13
    //
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt; ig++ )
        {
            s = - fe.phiDer( iloc, 2, ig ) * fe.phiDer( iloc, 0, ig ) * fe.weightDet( ig );
        }
        mat_tmp13( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 13
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s = - fe.phiDer( iloc, 2, ig ) * fe.phiDer( jloc, 0, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp13( iloc, jloc ) -= coef_s;
        mat_tmp13( jloc, iloc ) -= coef_s;
    }

    // diagonal 21
    //
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt; ig++ )
        {
            s = - fe.phiDer( iloc, 0, ig ) * fe.phiDer( iloc, 1, ig ) * fe.weightDet( ig );
        }
        mat_tmp21( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 21
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s = - fe.phiDer( iloc, 0, ig ) * fe.phiDer( jloc, 1, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp21( iloc, jloc ) -= coef_s;
        mat_tmp21( jloc, iloc ) -= coef_s;
    }

    // diagonal 22
    //
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt; ig++ )
        {
            s = fe.phiDer( iloc, 0, ig ) * fe.phiDer( iloc, 0, ig ) * fe.weightDet( ig )
                + fe.phiDer( iloc, 2, ig ) * fe.phiDer( iloc, 2, ig ) * fe.weightDet( ig ) ;
        }
        mat_tmp22( iloc, iloc ) += coef * s;
    }
    // extra diagonal 22
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
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
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt; ig++ )
        {
            s = - fe.phiDer( iloc, 2, ig ) * fe.phiDer( iloc, 1, ig ) * fe.weightDet( ig );
        }
        mat_tmp23( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 23
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s = - fe.phiDer( iloc, 2, ig ) * fe.phiDer( jloc, 1, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp23( iloc, jloc ) -= coef_s;
        mat_tmp23( jloc, iloc ) -= coef_s;
    }

    // diagonal 31
    //
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt; ig++ )
        {
            s = - fe.phiDer( iloc, 0, ig ) * fe.phiDer( iloc, 2, ig ) * fe.weightDet( ig );
        }
        mat_tmp31( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 31
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s = - fe.phiDer( iloc, 0, ig ) * fe.phiDer( jloc, 2, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp31( iloc, jloc ) -= coef_s;
        mat_tmp31( jloc, iloc ) -= coef_s;
    }

    // diagonal 32
    //
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt; ig++ )
        {
            s = - fe.phiDer( iloc, 1, ig ) * fe.phiDer( iloc, 2, ig ) * fe.weightDet( ig );
        }
        mat_tmp32( iloc, iloc ) -= coef * s;
    }
    // extra diagonal 32
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s = - fe.phiDer( iloc, 1, ig ) * fe.phiDer( jloc, 2, ig ) * fe.weightDet( ig );
        }
        coef_s = coef * s;
        mat_tmp32( iloc, jloc ) -= coef_s;
        mat_tmp32( jloc, iloc ) -= coef_s;
    }

    // diagonal 33
    //
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0; ig < fe.nbQuadPt; ig++ )
        {
            s = fe.phiDer( iloc, 0, ig ) * fe.phiDer( iloc, 0, ig ) * fe.weightDet( ig )
                + fe.phiDer( iloc, 1, ig ) * fe.phiDer( iloc, 1, ig ) * fe.weightDet( ig ) ;
        }
        mat_tmp33( iloc, iloc ) += coef * s;
    }
    // extra diagonal 33
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s = fe.phiDer( iloc, 1, ig ) * fe.phiDer( jloc, 1, ig ) * fe.weightDet( ig )
                + fe.phiDer( iloc, 2, ig ) * fe.phiDer( jloc, 2, ig ) * fe.weightDet( ig )       ;
        }
        coef_s = coef * s;
        mat_tmp33( iloc, jloc ) += coef_s;
        mat_tmp33( jloc, iloc ) += coef_s;
    }

    ElemMat::matrix_view mat_icomp = elmat.block( iblock + 0, jblock + 0 );
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) += mat_tmp11( iloc, iloc );
    }
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) += mat_tmp11( iloc, jloc );
        mat_icomp( jloc, iloc ) += mat_tmp11( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 0, jblock + 1 );
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp12( iloc, iloc );
    }
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp12( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp12( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 0, jblock + 2 );
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp13( iloc, iloc );
    }
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp13( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp13( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 1, jblock + 0 );
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp21( iloc, iloc );
    }
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp21( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp21( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 1, jblock + 1 );
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) += mat_tmp22( iloc, iloc );
    }
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) += mat_tmp22( iloc, jloc );
        mat_icomp( jloc, iloc ) += mat_tmp22( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 1, jblock + 2 );
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp23( iloc, iloc );
    }
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp23( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp23( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 2, jblock + 0 );
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp31( iloc, iloc );
    }
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp31( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp31( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 2, jblock + 1 );
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) -= mat_tmp32( iloc, iloc );
    }
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        mat_icomp( iloc, jloc ) -= mat_tmp32( iloc, jloc );
        mat_icomp( jloc, iloc ) -= mat_tmp32( jloc, iloc );
    }

    mat_icomp = elmat.block( iblock + 2, jblock + 2 );
    for ( i = 0; i < fe.nbDiag; i++ )
    {
        iloc = fe.patternFirst( i );
        mat_icomp( iloc, iloc ) += mat_tmp33( iloc, iloc );
    }
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
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

    ASSERT_PRE( fe.hasFirstDeriv(),
                "Stiffness div matrix needs at least the first derivatives" );
    double s;

    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( int icoor = 0; icoor < fe.nbCoor; ++icoor )
    {
        for ( int jcoor = 0; jcoor < fe.nbCoor; ++jcoor )
        {

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

            for ( int i = 0; i < fe.nbNode; ++i )
            {
                for ( int j = 0; j < fe.nbNode; ++j )
                {
                    s = 0;
                    for ( int ig = 0; ig < fe.nbQuadPt; ++ig )
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

    ASSERT_PRE( fe.hasFirstDeriv(),
                "Stiffness dergradbis matrix needs at least the first derivatives" );

    double s;
    Real guk[ fe.nbCoor ][ fe.nbCoor ][ fe.nbQuadPt ];      // \grad u^k at each quadrature point


    // loop on quadrature points
    for ( int ig = 0;ig < fe.nbQuadPt;ig++ )
    {

        // loop on space coordinates
        for ( int icoor = 0;icoor < fe.nbCoor;icoor++ )
        {

            // loop  on space coordinates
            for ( int jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            {
                s = 0.0;
                for ( int i = 0;i < fe.nbNode;i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbNode ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( int icoor = 0; icoor < fe.nbCoor; ++icoor )
    {
        for ( int jcoor = 0; jcoor < fe.nbCoor; ++jcoor )
        {

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

            for ( int i = 0; i < fe.nbNode; ++i )
            {
                for ( int j = 0; j < fe.nbNode; ++j )
                {
                    s = 0;
                    for ( int k = 0; k < fe.nbCoor; ++k )
                        for ( int ig = 0;ig < fe.nbQuadPt; ++ig )
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

    ASSERT_PRE( fe.hasFirstDeriv(),
                "Stiffness dergrad matrix needs at least the first derivatives" );

    double s;
    Real guk[ fe.nbCoor ][ fe.nbCoor ][ fe.nbQuadPt ];      // \grad u^k at each quadrature point


    // loop on quadrature points
    for ( int ig = 0;ig < fe.nbQuadPt;ig++ )
    {

        // loop on space coordinates
        for ( int icoor = 0;icoor < fe.nbCoor;icoor++ )
        {

            // loop  on space coordinates
            for ( int jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            {
                s = 0.0;
                for ( int i = 0;i < fe.nbNode;i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbNode ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //

    for ( int icoor = 0; icoor < fe.nbCoor; ++icoor )
    {
        for ( int jcoor = 0; jcoor < fe.nbCoor; ++jcoor )
        {

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

            for ( int i = 0; i < fe.nbNode; ++i )
            {
                for ( int j = 0; j < fe.nbNode; ++j )
                {
                    s = 0;
                    for ( int k = 0; k < fe.nbCoor; ++k )
                        for ( int ig = 0;ig < fe.nbQuadPt; ++ig )
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

    ASSERT_PRE( fe.hasFirstDeriv(),
                "stiff_derdiv needs at least the first derivatives" );

    Real guk[ fe.nbCoor ][ fe.nbCoor ][ fe.nbQuadPt ];      // \grad u^k at each quadrature point
    Real s;

    // loop on quadrature points
    for ( int ig = 0;ig < fe.nbQuadPt;ig++ )
    {

        // loop on space coordinates
        for ( int icoor = 0;icoor < fe.nbCoor;icoor++ )
        {

            // loop  on space coordinates
            for ( int jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            {
                s = 0.0;
                for ( int i = 0;i < fe.nbNode;i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbNode ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ][ ig ] = s;
            }
        }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( int icoor = 0; icoor < fe.nbCoor; ++icoor )
    {
        for ( int jcoor = 0; jcoor < fe.nbCoor; ++jcoor )
        {

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

            for ( int i = 0; i < fe.nbNode; ++i )
            {
                for ( int j = 0; j < fe.nbNode; ++j )
                {
                    s = 0;
                    for ( int k = 0; k < fe.nbCoor; ++k )
                        for ( int ig = 0;ig < fe.nbQuadPt;ig++ )
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
    ASSERT_PRE( fe.hasFirstDeriv(),
                "Stiffness Strain matrix needs at least the first derivatives" );
    double s;
    double tmp = coef * 0.5;

    ElemMat::matrix_type mat_tmp( fe.nbNode, fe.nbNode );

    for ( int i = 0; i < fe.nbNode; ++i )
    {
        for ( int j = 0; j < fe.nbNode; ++j )
        {
            s = 0;
            for ( int ig = 0; ig < fe.nbQuadPt; ++ig )
                for ( int icoor = 0; icoor < fe.nbCoor; ++icoor )
                    s += fe.phiDer( i, icoor, ig ) * fe.phiDer( j, icoor, ig ) * fe.weightDet( ig );
            mat_tmp( i, j ) = tmp * s;
        }
    }
    for ( int icoor = 0; icoor < fe.nbCoor; ++icoor )
    {
        ElemMat::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
    }

    for ( int icoor = 0; icoor < fe.nbCoor; ++icoor )
    {
        for ( int jcoor = 0; jcoor < fe.nbCoor; ++jcoor )
        {
            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );
            for ( int i = 0; i < fe.nbNode; ++i )
            {
                for ( int j = 0; j < fe.nbNode; ++j )
                {
                    s = 0;
                    for ( int ig = 0; ig < fe.nbQuadPt; ++ig )
                        s += fe.phiDer( i, jcoor, ig ) * fe.phiDer( j, icoor, ig ) * fe.weightDet( ig );
                    mat( i, j ) += tmp * s;
                }
            }
        }
    }
}



void mass_divw( Real coef, const ElemVec& w_loc, ElemMat& elmat, const CurrentFE& fe,
                int iblock, int jblock, int nb )
/*
  modified mass matrix: ( div w u,v )
*/
{

    ASSERT_PRE( fe.hasFirstDeriv(),
                "Mass matrix, (div w u, v) needs at least the first derivatives" );
    Tab2d mat_tmp( fe.nbNode, fe.nbNode );
    mat_tmp = ZeroMatrix( fe.nbNode, fe.nbNode );

    int i, icomp, ig, icoor, iloc, jloc;
    Real s, coef_s, divw[ fe.nbQuadPt ];

    // divw at quadrature nodes
    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {
        divw[ ig ] = 0.0;
        for ( icoor = 0; icoor < fe.nbCoor; ++icoor )
            for ( i = 0; i < fe.nbNode; ++i )
                divw[ ig ] += fe.phiDer( i, icoor, ig ) * w_loc.vec() [ i + icoor * fe.nbNode ];
    }

    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s += divw[ ig ] * fe.phi( iloc, ig ) * fe.phi( iloc, ig ) * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
            s += divw[ ig ] * fe.phi( iloc, ig ) * fe.phi( jloc, ig ) * fe.weightDet( ig );
        coef_s = coef * s;
        mat_tmp( iloc, jloc ) += coef_s;
        mat_tmp( jloc, iloc ) += coef_s;
    }
    // copy on the components
    for ( icomp = 0;icomp < nb;icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0;i < fe.nbDiag;i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
        {
            iloc = fe.patternFirst( i );
            jloc = fe.patternSecond( i );
            mat_icomp( iloc, jloc ) += mat_tmp( iloc, jloc );
            mat_icomp( jloc, iloc ) += mat_tmp( jloc, iloc );
        }
    }
}


void mass_divw(const std::vector<Real>& coef, const ElemVec& w_loc, ElemMat& elmat, const CurrentFE& fe,
                int iblock, int jblock, int nb )
/*
  modified mass matrix: ( div w u,v )
*/
{

    ASSERT_PRE( fe.hasFirstDeriv(),
                "Mass matrix, (div w u, v) needs at least the first derivatives" );
    Tab2d mat_tmp( fe.nbNode, fe.nbNode );
    mat_tmp = ZeroMatrix( fe.nbNode, fe.nbNode );

    int i, icomp, ig, icoor, iloc, jloc;
    Real s, divw[ fe.nbQuadPt ]; // , coef_s

    // divw at quadrature nodes
    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {
        divw[ ig ] = 0.0;
        for ( icoor = 0; icoor < fe.nbCoor; ++icoor )
            for ( i = 0; i < fe.nbNode; ++i )
                divw[ ig ] += fe.phiDer( i, icoor, ig ) * w_loc.vec() [ i + icoor * fe.nbNode ];
    }

    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            s += coef[ig] * divw[ ig ] * fe.phi( iloc, ig ) * fe.phi( iloc, ig ) * fe.weightDet( ig );
        }
        mat_tmp( iloc, iloc ) += s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
            s += coef[ig]* divw[ ig ] * fe.phi( iloc, ig ) * fe.phi( jloc, ig ) * fe.weightDet( ig );
        mat_tmp( iloc, jloc ) += s;
        mat_tmp( jloc, iloc ) += s;
    }
    // copy on the components
    for ( icomp = 0;icomp < nb;icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0;i < fe.nbDiag;i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
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

    ASSERT_PRE( fe.hasFirstDeriv(),
                "Mass matrix, (grad u_0 u, v) needs at least the first derivatives" );

    int ig, icoor, jcoor, i, j;
    Real s;
    Real gu0[ fe.nbQuadPt ][ fe.nbCoor ][ fe.nbCoor ];


    //
    // grad u0 at quadrature nodes
    //
    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {
        for ( icoor = 0; icoor < fe.nbCoor; ++icoor )
            for ( jcoor = 0; jcoor < fe.nbCoor; ++jcoor )
            {
                gu0[ ig ][ icoor ][ jcoor ] = 0.0;
                for ( i = 0; i < fe.nbNode; ++i )
                    gu0[ ig ][ icoor ][ jcoor ] += fe.phiDer( i, jcoor, ig ) * u0_loc.vec() [ i + icoor * fe.nbNode ];
            }
    }
    //
    // blocks (icoor,jcoor) of elmat
    //
    for ( icoor = 0; icoor < fe.nbCoor; ++icoor )
    {
        for ( jcoor = 0; jcoor < fe.nbCoor; ++jcoor )
        {

            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );

            for ( i = 0; i < fe.nbNode; ++i )
            {
                for ( j = 0; j < fe.nbNode; ++j )
                {
                    s = 0;
                    for ( ig = 0; ig < fe.nbQuadPt; ++ig )
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
    ASSERT_PRE( fe.hasFirstDeriv(),
                "Stiffness matrix needs at least the first derivatives" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int iloc, jloc;
    int i, icoor, ig, jcoor;
    double s, coef_s, coef_v[ nDimensions ];
    //    int nbN1=fe.nbNode;
    int nbN2 = fe2.nbNode;
    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            for ( icoor = 0;icoor < (int)nDimensions;icoor++ )
                coef_v[ icoor ] = 0.;

            // computation of the convection term in the quadrature nodes
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
            {
                for ( int iloc = 0;iloc < nbN2;iloc++ )
                    coef_v[ icoor ] += vec_loc.vec() [ iloc + icoor * nbN2 ] * fe2.phi( iloc, ig );
            }

            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
            {
                for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
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
    for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            for ( icoor = 0;icoor < (int)nDimensions;icoor++ )
                coef_v[ icoor ] = 0.;

            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
            {
                for ( int iloc = 0;iloc < nbN2;iloc++ )
                    coef_v[ icoor ] += vec_loc.vec() [ iloc + icoor * nbN2 ] * fe2.phi( iloc, ig );
            }

            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
            {
                for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
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
    for ( int icomp = 1;icomp < nb;icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0;i < fe.nbDiag;i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat( iloc, iloc );
        }
        for ( i = fe.nbDiag;i < fe.nbDiag + fe.nbUpper;i++ )
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
    ASSERT_PRE( fe_u.hasFirstDeriv(),
                "Gradient matrix needs at least the first derivatives" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int ig;
    int i, j;
    double s;
    for ( i = 0;i < fe_u.nbNode;i++ )
    {
        for ( j = 0;j < fe_p.nbNode;j++ )
        {
            s = 0;
            for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )
                // Be careful the minus is here and not in coef!!!!

                // wrong for different quadrules of fe_u and fe_p !!   Martin P.
                // s -= fe_p.phi(j,ig)*fe_u.phiDer(i,icoor,ig)*fe_u.weightDet(ig);

                s -= fe_p.refFE.phi( j, fe_u.qr.quadPointCoor( ig, 0 ), fe_u.qr.quadPointCoor( ig, 1 ),
                                     fe_u.qr.quadPointCoor( ig, 2 ) ) * fe_u.phiDer( i, icoor, ig ) * fe_u.weightDet( ig );
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
    ASSERT_PRE( fe_u.hasFirstDeriv(),
                "Divergence matrix needs at least the first derivatives" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int ig;
    int i, j;
    double s;
    for (i = 0; i < fe_u.nbNode; i++)
    {
        for (j = 0; j < fe_p.nbNode; j++)
        {
            s = 0;
            for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )
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
    ASSERT_PRE( fe_u.hasFirstDeriv(),
                "grad_div: Gradient matrix needs at least the first derivatives" );
    double s;
    int iblock = block_pres - nDimensions;
    for ( int icoor = 0;icoor < 3;icoor++ )
    {
        ElemMat::matrix_view mat_grad = elmat.block( iblock + icoor, block_pres );
        ElemMat::matrix_view mat_div = elmat.block( block_pres , iblock + icoor );
        for ( int i = 0;i < fe_u.nbNode;i++ )
        {
            for ( int j = 0;j < fe_p.nbNode;j++ )
            {
                s = 0;
                for ( int ig = 0;ig < fe_u.nbQuadPt;ig++ )
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
    for ( int i = 0;i < fe.nbNode;i++ )
    {
        for ( int j = 0;j < fe.nbNode;j++ )
        {
            s = 0;
            for ( int ig = 0;ig < fe.nbQuadPt;ig++ )
            {
                for ( int icoor = 0;icoor < 3;icoor++ )
                {
                    s += fe.phiDer( i, icoor, ig ) * fe.phiDer( j, icoor, ig ) * fe.weightDet( ig );
                }
            }
            mat( i, j ) -= fh2 * s;
        }
    }
}

void advection( Real /*coef*/, ElemVec& vel,
                ElemMat& elmat, const CurrentFE& fe, int iblock, int jblock, int nb )
{
    ASSERT_PRE( fe.hasFirstDeriv(),
                "advection (vect) matrix needs at least the first derivatives" );
    //ASSERT_PRE(nb>1,"if nb = 1, use the other advection function");

    Tab2d mat_tmp( fe.nbNode, fe.nbNode );
    Real v_grad, s;
    Real v[ nDimensions ];
    for ( int ig = 0;ig < fe.nbQuadPt;ig++ )
    {
        for ( int icoor = 0;icoor < ( int ) nDimensions;icoor++ )
        {
            ElemVec::vector_view velicoor = vel.block( icoor );
            v[ icoor ] = 0.;
            for ( int k = 0;k < fe.nbNode;k++ )
            {
                v[ icoor ] += velicoor( k ) * fe.phi( k, ig ); // velocity on the intgt point
            }
        }
        for ( int i = 0;i < fe.nbNode;i++ )
        {
            for ( int j = 0;j < fe.nbNode;j++ )
            {
                s = 0.;
                for ( int ig = 0;ig < fe.nbQuadPt;ig++ )
                {
                    v_grad = 0.;
                    for ( int icoor = 0;icoor < ( int ) nDimensions;icoor++ )
                    {
                        v_grad += v[ icoor ] * fe.phiDer( j, icoor, ig );
                    }
                    s += v_grad * fe.phi( i, ig ) * fe.weightDet( ig );
                }
                mat_tmp( i, j ) = s;
            }
        }
    }
    // copy on the components
    for ( int icomp = 0;icomp < nb;icomp++ )
    {
        ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( int i = 0;i < fe.nbDiag;i++ )
        {
            for ( int j = 0;j < fe.nbDiag;j++ )
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
        int iq;
        int i, j;
        double s, coef;
        int nbN1 = fe1.nbNode;
        for ( i = 0;i < nbN1;i++ )
        {
            for ( j = 0;j < fe2.nbNode;j++ )
            {
                s = 0;
                for ( iq = 0;iq < fe1.nbQuadPt;iq++ )
                {
                    coef = 0;
                    for ( int iloc = 0;iloc < nbN1;iloc++ )
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
        int iq;
        int i, j;
        double s, coef, coef_div;
        int nbN1 = fe1.nbNode;
        for ( i = 0;i < nbN1;i++ )
        {
            for ( j = 0;j < fe2.nbNode;j++ )
            {
                s = 0;
                for ( iq = 0;iq < fe1.nbQuadPt;iq++ )
                {
                    coef = 0;
                    coef_div = 0;

                    for ( int iloc = 0;iloc < nbN1;iloc++ )
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
        int iq;
        int i, j;
        double s, coef;
        int nbN1 = fe1.nbNode;
        int nbN3 = fe3.nbNode;

        for ( i = 0; i < nbN1; i++ )
        {
            for ( j = 0; j < fe2.nbNode; j++ )
            {
                s = 0;
                for ( iq = 0; iq < fe1.nbQuadPt; iq++ )
                {
                    coef = 0;

                    for ( int iloc = 0;iloc < nbN3; iloc++ )
                        coef += vec_loc.vec() [iloc + icoor*nbN3] * fe3.phi(iloc, iq);

                    s += coef * fe2.phi(i, iq) * fe1.phiDer(j, icoor, iq) * fe1.weightDet(iq);
                } // Loop on quadrature nodes

                mat(i, j) += s;
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
for(i=0;i<fe.nbDiag;i++){
iloc = fe.patternFirst(i);
s = 0;
for (iq=0;i<siz;++i){
qloc = fe.patternFirst(iq);
for(ig=0;ig<fe.nbQuadPt;ig++)
s += coef(iq)*fe.phi(qloc,ig)*fe.phi(iloc,ig)*fe.phi(iloc,ig)*
fe.weightDet(ig);
}
mat(iloc,iloc) += 2*s;
elvec(iloc) += coef(iloc)*s;
}
//
// extra diagonal
//
for(i=fe.nbDiag;i<fe.nbDiag+fe.nbUpper;i++){
iloc = fe.patternFirst(i);
jloc = fe.patternSecond(i);
s = 0;
for (iq=0;i<siz;++i){
qloc = fe.patternFirst(iq);
for(ig=0;ig<fe.nbQuadPt;ig++)
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
    int i, ig;
    ASSERT_PRE( fe.hasJac(), "Source vector needs at least the jacobian" );
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;
    for ( i = 0;i < fe.nbNode;i++ )
    {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
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
    int i, ig;
    ASSERT_PRE( fe.hasJac(), "Source vector (f) needs at least the jacobian" );
    ElemVec::vector_view vec = elvec.block( eblock );
    ElemVec::vector_view vecf = f.block( fblock );
    Real f_ig;

    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {
        f_ig = 0.;
        for ( i = 0;i < fe.nbNode;i++ )
            f_ig += vecf( i ) * fe.phi( i, ig );
        for ( i = 0;i < fe.nbNode;i++ )
        {
            vec( i ) += coef * f_ig * fe.phi( i, ig ) * fe.weightDet( ig );
        }
    }
}




  void source_divuq(Real alpha, ElemVec& uLoc,  ElemVec& elvec, const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock  )
  {
    int i, j, ic, iq;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;

    for (i = 0; i < fe_p.nbNode; i++) {
      s = 0;
      for (iq = 0; iq < fe_p.nbQuadPt; ++iq )
	for (j = 0; j < fe_u.nbNode; ++j)
	  for (ic = 0; ic < (int)nDimensions; ++ic)
	    s += uLoc[ic*fe_u.nbNode+j]*fe_u.phiDer(j,ic,iq)*fe_p.phi(i,iq)* fe_p.weightDet( iq );

      vec( i ) += s*alpha;
    }

  }

  void source_gradpv(Real alpha, ElemVec& pLoc,  ElemVec& elvec, const CurrentFE& fe_p, const CurrentFE& fe_u, int iblock )
  {
    int i, j, iq;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;

    for ( i = 0;i < fe_u.nbNode;i++ )
      {
        s = 0;
        for (iq = 0; iq < fe_u.nbQuadPt; ++iq )
	  for (j = 0; j < fe_p.nbNode; ++j)
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
    int i, ig;
    ASSERT_PRE( fe.hasJac(), "elemOper.cpp source_fhn: FHN source vector (f) needs at least the jacobian" );
    ElemVec::vector_view vec = elvec.block( eblock );
    ElemVec::vector_view vecu = u.block( fblock );
    Real f_ig;

    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {
        f_ig = 0.;
        for ( i = 0;i < fe.nbNode;i++ )
            f_ig += vecu( i ) * fe.phi( i, ig );
        for ( i = 0;i < fe.nbNode;i++ )
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
    ASSERT_PRE( fe.hasFirstDeriv(),
                "source_mass needs at least the first derivatives" );


    Real B[ fe.nbCoor ][ fe.nbCoor ];                 // \grad (convect) at a quadrature point
    Real A[ fe.nbCoor ][ fe.nbCoor ];                 // I\div d - (\grad d)^T at a quadrature point
    Real aux[ fe.nbQuadPt ];                        // grad (- w^k):[I\div d - (\grad d)^T] at  quadrature points
    Real uk[ fe.nbQuadPt ][ fe.nbCoor ];              // u^k quadrature points
    Real guk[ fe.nbQuadPt ][ fe.nbCoor ][ fe.nbCoor ];  // \grad u^k at quadrature points
    Real convect[ fe.nbCoor ];                      // convect at quadrature points
    Real convect_A[ fe.nbQuadPt ][ fe.nbCoor ];       // (convect)^T [I\div d - (\grad d)^T] at quadrature points

    Real s, sA, sB, sG;


    int icoor, jcoor, ig, i;
    // loop on quadrature points
    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {

        // loop on space coordindates
        for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
        {

            // each compontent of uk at each quadrature points
            s = 0.0;
            for ( i = 0;i < fe.nbNode;i++ )
                s += fe.phi( i, ig ) * uk_loc.vec() [ i + icoor * fe.nbNode ];
            uk[ ig ][ icoor ] = s;//uk_x(pt_ig), uk_y(pt_ig), uk_z(pt_ig)

            // each compontent of convect at this quadrature point
            s = 0.0;
            for ( i = 0;i < fe.nbNode;i++ )
                s += fe.phi( i, ig ) * convect_loc.vec() [ i + icoor * fe.nbNode ];
            convect[ icoor ] = s;


            // loop  on space coordindates
            for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            {
                sB = 0.0;
                sA = 0.0;
                sG = 0.0;
                for ( i = 0;i < fe.nbNode;i++ )
                {
                    sG += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbNode ]; //  \grad u^k at each quadrature point
                    sB -= fe.phiDer( i, jcoor, ig ) * wk_loc.vec() [ i + icoor * fe.nbNode ]; //  \grad (- w^k) at this quadrature point
                    sA -= fe.phiDer( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe.nbNode ]; //  - (\grad d) ^T at this quadrature point
                }
                guk[ ig ][ icoor ][ jcoor ] = sG; // \grad u^k at each quadrature point
                B[ icoor ][ jcoor ] = sB; // \grad (convect) at this quadrature point
                A[ icoor ][ jcoor ] = sA; // -(\grad d) ^T at this quadrature point
            }
        }

        s = 0.0;
        for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            s -= A[ jcoor ][ jcoor ];  // \div d at this quadrature point ( - trace( A ) )

        for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            A[ jcoor ][ jcoor ] += s;  // I\div d - (\grad d)^T at this quadrature point (A+I(-tr(A)))

        s = 0;
        for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
            for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
                s += B[ icoor ][ jcoor ] * A[ icoor ][ jcoor ]; // \grad (-w^k):[I\div d - (\grad d)^T] at each quadrature point
        aux[ ig ] = s;

        s = 0;
        for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
        {
            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                s += convect[ icoor ] * A[ icoor ][ jcoor ]; // convect^T [I\div d - (\grad d)^T]
            convect_A[ ig ][ jcoor ] = s;
        }
//         std::cout<<"aux "<<aux[ig]<<std::endl;

//     for(icoor=0; icoor<fe.nbCoor; ++icoor)
//         for(jcoor=0; jcoor<fe.nbCoor; ++jcoor)
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
    for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
    {

        // the block iccor of the elementary vector
        ElemVec::vector_view vec = elvec.block( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0;i < fe.nbNode;i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0;ig < fe.nbQuadPt;ig++ )
            {

                // \grad ( - w^k ):[I\div d - (\grad d)^T] \phi_i
                s += aux[ ig ] * uk[ ig ][ icoor ] * fe.phi( i, ig ) * fe.weightDet( ig );

                // convect^T [I\div d - (\grad d)^T] (\grad u^k)^T \phi_i
                for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
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

    ASSERT_PRE( fe.hasFirstDeriv(),
                "source_mass needs at least the first derivatives" );

    Real guk[ fe.nbCoor ][ fe.nbCoor ];      // \grad u^k at a quadrature point
    Real dw[ fe.nbCoor ];                  // dw at a quadrature point
    Real aux[ fe.nbQuadPt ][ fe.nbCoor ];    // (\grad u^k)dw at each quadrature point
    Real s;

    int ig, icoor, jcoor, i;

    // loop on quadrature points
    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
        {

            // each compontent (icoor) of dw at this quadrature point
            s = 0.0;
            for ( i = 0;i < fe.nbNode;i++ )
                s += fe.phi( i, ig ) * dw_loc.vec() [ i + icoor * fe.nbNode ];
            dw[ icoor ] = s;

            // loop  on space coordinates
            for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            {
                s = 0.0;
                for ( i = 0;i < fe.nbNode;i++ )
                    s += fe.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe.nbNode ]; //  \grad u^k at a quadrature point
                guk[ icoor ][ jcoor ] = s;
            }
        }

        // (\grad u^k)dw at each quadrature point
        for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
        {
            s = 0.0;
            for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
                s += guk[ icoor ][ jcoor ] * dw[ jcoor ];
            aux[ ig ][ icoor ] = s;
        }
    }

    //
    // Numerical integration
    //

    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
    {

        ElemVec::vector_view vec = elvec.block( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0;i < fe.nbNode;i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0;ig < fe.nbQuadPt;ig++ )
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
    ASSERT_PRE( fe.hasFirstDeriv(),
                "source_mass needs at least the first derivatives" );


    Real B[ fe.nbCoor ][ fe.nbCoor ];                 // \grad u^n at a quadrature point
    Real A[ fe.nbCoor ][ fe.nbCoor ];                 // I\div d - (\grad d)^T at a quadrature point
    Real aux[ fe.nbQuadPt ];                        //  \div d  \div u^n  + grad u^n:[I\div d - (\grad d)^T] at  quadrature points
    Real uk[ fe.nbQuadPt ][ fe.nbCoor ];              // u^k quadrature points

    Real s, sA, sB;


    int icoor, jcoor, ig, i;
    // loop on quadrature points
    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {

        // loop on space coordindates
        for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
        {

            // each compontent of uk at each quadrature points
            s = 0.0;
            for ( i = 0;i < fe.nbNode;i++ )
                s += fe.phi( i, ig ) * uk_loc.vec() [ i + icoor * fe.nbNode ];
            uk[ ig ][ icoor ] = s;


            // loop  on space coordindates
            for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            {
                sB = 0.0;
                sA = 0.0;
                for ( i = 0;i < fe.nbNode;i++ )
                {
                    sB += fe.phiDer( i, jcoor, ig ) * un_loc.vec() [ i + icoor * fe.nbNode ]; //  \grad u^n at this quadrature point
                    sA -= fe.phiDer( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe.nbNode ]; //  - (\grad d) ^T at this quadrature point
                }
		B[ icoor ][ jcoor ] = sB; // \grad u^n at this quadrature point
                A[ icoor ][ jcoor ] = sA; // -(\grad d) ^T at this quadrature point
            }
        }

        s = 0.0;
        for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
	  s -= A[ jcoor ][ jcoor ];  // \div d at this quadrature point ( - trace( A ) )

        for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
	  A[ jcoor ][ jcoor ] += 2 * s;  // 2 * I\div d - (\grad d)^T at this quadrature point

        s = 0;
        for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
            for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
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
    for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
    {

      // the block iccor of the elementary vector
      ElemVec::vector_view vec = elvec.block( icoor );

      // loop on nodes, i.e. loop on components of this block
      for ( i = 0;i < fe.nbNode;i++ )
        {
	  // loop on quadrature points
	  s = 0;
	  for ( ig = 0;ig < fe.nbQuadPt;ig++ )
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

    ASSERT_PRE( fe_u.hasFirstDeriv(),
                "source_stress needs at least the velocity shape functions first derivatives" );

    Real A[ fe_u.nbCoor ][ fe_u.nbCoor ];                 // I\div d - (\grad d)^T at a quadrature point
    Real guk[ fe_u.nbCoor ][ fe_u.nbCoor ];               // \grad u^k at a quadrature point
    Real sigma[ fe_u.nbCoor ][ fe_u.nbCoor ];             // [-p^k I + 2*mu e(u^k)] a quadrature point
    Real B[ fe_u.nbQuadPt ][ fe_u.nbCoor ][ fe_u.nbCoor ];  // [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] at each quadrature point
    Real s, sA, sG, pk;

    int icoor, jcoor, kcoor, ig, i;

    // loop on quadrature points
    for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
        {

            // loop  on space coordindates
            for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
            {
                sA = 0.0;
                sG = 0.0;
                for ( i = 0;i < fe_u.nbNode;i++ )
                {
                    sG += fe_u.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbNode ]; //  \grad u^k at this quadrature point
                    sA -= fe_u.phiDer( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe_u.nbNode ]; //  - (\grad d) ^T at this quadrature point
                }
                guk[ icoor ][ jcoor ] = sG;
                A[ icoor ][ jcoor ] = sA;
            }
        }

        s = 0.0;
        for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
            s -= A[ jcoor ][ jcoor ];  // \div d at a quadrature point ( - trace( A ) )

//         std::cout<<"div = "<< s <<std::endl;

        for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
            A[ jcoor ][ jcoor ] += s;  // I\div d  - (\grad d)^T

        pk = 0.0;
        for ( i = 0;i < fe_p.nbNode;i++ )
            pk += fe_p.phi( i, ig ) * pk_loc.vec() [ i ]; // p^k at this quadrature point


        // sigma = [-p^k I + 2*mu e(u^k)] a quadrature point
        for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
        {
            for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
                sigma[ icoor ][ jcoor ] = mu * ( guk[ icoor ][ jcoor ] + guk[ jcoor ][ icoor ] );
            sigma[ icoor ][ icoor ] -= pk;
        }

        // [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] at each quadrature point
        for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
            for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
            {
                s = 0;
                for ( kcoor = 0;kcoor < fe_u.nbCoor;kcoor++ )
                    s += sigma[ icoor ][ kcoor ] * A[ kcoor ][ jcoor ];
                B[ ig ][ icoor ][ jcoor ] = s;
            }
     }

//     for(icoor=0; icoor<fe_u.nbCoor; ++icoor)
//         for(jcoor=0; jcoor<fe_u.nbCoor; ++jcoor)
//             {
//                 std::cout<<" Atrue ["<<icoor<<"]["<<jcoor<<"] = "<<A[icoor][jcoor]<<std::endl;
//             }

//     for(icoor=0; icoor<fe_u.nbCoor; ++icoor)
//         for(jcoor=0; jcoor<fe_u.nbCoor; ++jcoor)
//             {
//                 double l=0.;
//                 for(int e=0; e<fe_u.nbQuadPt; ++e)
//                     l+=B[e][icoor][jcoor];
//                 std::cout<<" Btrue ["<<icoor<<"]["<<jcoor<<"] = "<<l<<std::endl;
//             }

    //
    // Numerical integration
    //

    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
    {

        ElemVec::vector_view vec = elvec.block( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0;i < fe_u.nbNode;i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )

                for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
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

    ASSERT_PRE( fe_u.hasFirstDeriv(),
                "source_stress needs at least the velocity shape functions first derivatives" );


    Real guk[ fe_u.nbCoor ][ fe_u.nbCoor ];               // \grad u^k at a quadrature point
    Real gd[ fe_u.nbCoor ][ fe_u.nbCoor ];                // \grad d at a quadrature point
    Real A[ fe_u.nbQuadPt ][ fe_u.nbCoor ][ fe_u.nbCoor ];  // \grad u^k \grad d + [\grad d]^T[\grad u^k]^T  at each quadrature point
    Real su, sd, s;

    int icoor, jcoor, kcoor, ig, i;

    // loop on quadrature points
    for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
        {

            // loop  on space coordindates
            for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
            {
                su = 0.0;
                sd = 0.0;
                for ( i = 0;i < fe_u.nbNode;i++ )
                {
                    su += fe_u.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbNode ]; //  \grad u^k at this quadrature point
                    sd += fe_u.phiDer( i, jcoor, ig ) * d_loc.vec() [ i + icoor * fe_u.nbNode ];  //  \grad d at this quadrature point
                }
                guk[ icoor ][ jcoor ] = su;
                gd[ icoor ][ jcoor ] = sd;
            }
        }


        for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
        {
            for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
            {
                s = 0;
                for ( kcoor = 0;kcoor < fe_u.nbCoor;kcoor++ )
                    s += guk[ icoor ][ kcoor ] * gd[ kcoor ][ jcoor ] + gd[ kcoor ][ icoor ] * guk[ jcoor ][ kcoor ];
                A[ ig ][ icoor ][ jcoor ] = s;
            }
        }

//     for(icoor=0; icoor<fe_u.nbCoor; ++icoor)
//         for(jcoor=0; jcoor<fe_u.nbCoor; ++jcoor)
//             {
//                 std::cout<<" gdtrue ["<<icoor<<"]["<<jcoor<<"] = "<<gd[icoor][jcoor]<<std::endl;
//             }
//     for(icoor=0; icoor<fe_u.nbCoor; ++icoor)
//         for(jcoor=0; jcoor<fe_u.nbCoor; ++jcoor)
//             {
//                 std::cout<<" Atrue ["<<icoor<<"]["<<jcoor<<"] = "<<A[ig][icoor][jcoor]<<std::endl;
//             }
    }

    //
    // Numerical integration
    //
    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
    {

        ElemVec::vector_view vec = elvec.block( icoor );

        // loop on nodes, i.e. loop on components of this block
        for ( i = 0;i < fe_u.nbNode;i++ )
        {

            // loop on quadrature points
            s = 0;
            for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )

                for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
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

    ASSERT_PRE( fe_u.hasFirstDeriv(),
                "source_stress needs at least the velocity shape functions first derivatives" );
    Real A[ fe_u.nbCoor ][ fe_u.nbCoor ];     //  I\div d - (\grad d)^T at a quadrature point
    Real guk[ fe_u.nbCoor ][ fe_u.nbCoor ];   // \grad u^k at a quadrature point
    Real aux[ fe_u.nbQuadPt ];              // grad u^k:[I\div d - (\grad d)^T] at each quadrature point
    ElemVec::vector_view vec = elvec.block( iblock );

    Real s, sA, sG;
    int icoor, jcoor, ig, i;


    // loop on quadrature points
    for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )
    {

        // loop on space coordinates
        for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
        {

            // loop  on space coordinates
            for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
            {
                sA = 0.0;
                sG = 0.0;
                for ( i = 0;i < fe_u.nbNode;i++ )
                {
                    sG += fe_u.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbNode ]; //  \grad u^k at a quadrature point
                    sA -= fe_u.phiDer( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe_u.nbNode ]; //  - (\grad d) ^T at a quadrature point
                }
                guk[ icoor ][ jcoor ] = sG;
                A[ icoor ][ jcoor ] = sA;
            }
        }

        s = 0.0;
        for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
            s -= A[ jcoor ][ jcoor ];  // \div d at this quadrature point ( - trace( A ) )

        for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
            A[ jcoor ][ jcoor ] += s;  // I\div d - (\grad d)^T at this quadrature point

        s = 0;
        for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
            for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
                s += guk[ icoor ][ jcoor ] * A[ icoor ][ jcoor ]; // \grad u^k : [I\div d - (\grad d)^T] at each quadrature point
        aux[ ig ] = s;
    }

    //
    // Numerical integration
    //

    // Loop on nodes, i.e. loop on elementary vector components
    for ( i = 0;i < fe_p.nbNode;i++ )
    {

        // loop on quadrature points
        s = 0;
        for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )
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

  ASSERT_PRE( fe.hasFirstDeriv(),
	      "source_stress needs at least the velocity shape functions first derivatives" );

  Real A[ fe.nbCoor ][ fe.nbCoor ];     //  I\div d - (\grad d)^T - \grad d at a quadrature point
  Real B[ fe.nbCoor ][ fe.nbCoor ];     // - \grad d
  Real gpk[ fe.nbCoor ];   // \grad p^k at a quadrature point
  Real aux[ fe.nbQuadPt ][fe.nbCoor];              // [I\div d - (\grad d)^T - \grad d ]\grad p^k at each quadrature point

  ElemVec::vector_view vec = elvec.block( iblock );

  Real s, sA, sG;
  int icoor, jcoor, ig, i;


  // loop on quadrature points
  for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {


      // loop on space coordinates
      for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
        {

	  sG = 0.0;
	  for ( i = 0;i < fe.nbNode;i++ )
	    sG += fe.phiDer( i, icoor, ig ) * p_loc.vec() [ i ]; //  \grad p^k at a quadrature point
	  gpk[ icoor ] = sG;

	  // loop  on space coordinates
	  for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            {
	      sA = 0.0;
	      for ( i = 0;i < fe.nbNode;i++ )
		sA -= fe.phiDer( i, icoor, ig ) * d_loc.vec() [ i + jcoor * fe.nbNode ]; //  - (\grad d) ^T at a quadrature point
	      A[ icoor ][ jcoor ] = sA;
	      B[ jcoor ][ icoor ] = sA;
            }
        }

        s = 0.0;
        for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            s -= A[ jcoor ][ jcoor ];  // \div d at this quadrature point ( - trace( A ) )

        for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            A[ jcoor ][ jcoor ] += s;  // I\div d - (\grad d)^T at this quadrature point

	for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
	  for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
	    A[ icoor ][ jcoor ] += B[ icoor ][ jcoor ]; // I\div d - (\grad d)^T - \grad d at this quadrature point

        s = 0;
        for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
	  {
            for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
	      s += A[ icoor ][ jcoor ] * gpk[jcoor]; // [I\div d - (\grad d)^T -\grad d] \grad p^k at each quadrature point
	    aux[ ig ][icoor] = s;
	  }
    }

    //
    // Numerical integration
    //

    // Loop on nodes, i.e. loop on elementary vector components
    for ( i = 0;i < fe.nbNode;i++ )
    {

        // loop on quadrature points
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
	  for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
            s += aux[ ig ][icoor] * fe.phiDer( i, icoor, ig ) * fe.weightDet( ig );
        vec [ i ] += coef * s;
    }
}



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


//
//Cholesky decomposition
void choldc( KNM<Real> &a, KN<Real> &p )
{
    int i, j, k;
    Real sum;

    int n = a.N();
    for ( i = 0;i < n;i++ )
    {
        for ( j = i;j < n;j++ )
        {
            for ( sum = a( i, j ), k = i - 1;k >= 0;k-- )
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
    for ( i = 0;i < n;i++ )
    {
        for ( sum = b( i ), k = i - 1;k >= 0;k-- )
            sum -= a( i, k ) * x( k );
        x( i ) = sum / p( i );
    }
    for ( i = n - 1;i >= 0;i-- )
    {
        for ( sum = x( i ), k = i + 1;k < n;k++ )
            sum -= a( k, i ) * x( k );
        x( i ) = sum / p( i );
    }
}


void source_press( Real coef, const ElemVec& uk_loc, ElemMat& elmat,
                   const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock )
{

    ASSERT_PRE( fe_u.hasFirstDeriv(),
                "source_stress needs at least the velocity shape functions first derivatives" );
    Real A[ fe_u.nbCoor ][ fe_u.nbCoor ][fe_u.nbNode][ fe_u.nbCoor ];     //  I\div d - (\grad d)^T at a quadrature point
    Real guk[ fe_u.nbCoor ][ fe_u.nbCoor ];   // \grad u^k at a quadrature point
    Real aux[ fe_u.nbQuadPt ][fe_u.nbNode][ fe_u.nbCoor ];              // grad u^k:[I\div d - (\grad d)^T] at each quadrature point

    //    double A[ fe_u.nbCoor ][ fe_u.nbCoor ][fe_u.nbNode][ fe_u.nbCoor ];

    Real /*s,*/ l/*[fe_u.nbNode][fe_u.nbCoor]*/, sA/*[fe_u.nbNode][fe_u.nbCoor]*/, sG;
    int icoor, jcoor, ig;


    // loop on quadrature points
    for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )
        {


            ////INIT
            for(int p=0; p<fe_u.nbCoor; ++p)
                {
                    for(int q=0; q<fe_u.nbCoor; ++q)
                        {
                            for(int d=0; d<fe_u.nbNode; ++d)
                                {
                                    for(int e=0; e<fe_u.nbCoor; ++e)
                                        A[p][q][d][e]=0.;
                                }
                            //guk[p][q]=0.;
                        }
                    //for(int h=0; h<fe_u.nbNode; ++h)
                        //aux[ig][h][p]=0.;
                }
    ////END INIT

        // loop on space coordinates
        for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
        {

            //for ( i = 0;i < fe_u.nbNode;i++ )
                //for ( short k = 0;k < fe_u.nbCoor;k++ )
                    //sA[i][k] = 0.0;
            // loop  on space coordinates
            for ( UInt kcoor = 0;(int)kcoor < fe_u.nbCoor;kcoor++ )
            for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
            {
                sG = 0.0;
                for ( int i = 0;i < fe_u.nbNode;i++ )
                {
                    sG += fe_u.phiDer( i, jcoor, ig ) * uk_loc.vec() [ i + icoor * fe_u.nbNode ]; //  \grad u^k at a quadrature point
                    sA = -fe_u.phiDer( i, icoor, ig )/**fe_u.phi( i, ig )*/ ;/** d_loc.vec() [ i + jcoor * fe_u.nbNode ];*/ //  - (\grad d) ^T at a quadrature point
                    A[ icoor ][ jcoor ][i][jcoor] = sA; // -\delta_{jcoor kcoor} \partial_{icoor}
                }
                guk[ icoor ][ jcoor ] = sG;
            }
        }





        double z[fe_u.nbCoor];
        for ( int i = 0;i < fe_u.nbNode;i++ )
            {
                //                for ( int kcoor = 0;kcoor < fe_u.nbCoor;kcoor++ )
                for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
                    {
                                //if(icoor==jcoor)
                                z[ jcoor ] = A[ jcoor ][ jcoor ][i][jcoor];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
                    }

                for ( jcoor = 0;jcoor < fe_p.nbCoor;jcoor++ )
                    {
                        for (int kcoor = 0;kcoor < fe_p.nbCoor;kcoor++ )
                        {
                                //if(icoor==jcoor)
                                A[ jcoor ][ jcoor ][i][kcoor] -= z[ kcoor ];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
                        }
                    }



//                         for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
//                             {
//                             //if(kcoor==jcoor)
//                                     for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
//                                         {
//                                             //    if(icoor==jcoor)
//                                             A[ icoor ][ jcoor ][i][kcoor] -= A[ kcoor ][ icoor ][i][jcoor];  // \delta_{jcoor icoor}\partial_{kcoor}
//                                             //else
//                                             //f[icoor][jcoor]=0.;
//                                         }
//                             }
                        //                    }
                        //                for ( UInt kcoor = 0;kcoor < fe_u.nbCoor;kcoor++ )
                        //                    {
//                         for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
//                         {
//                             for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
//                                 {
//                                     //for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
//                                     if(icoor==jcoor)
//                                         A[ icoor ][ jcoor ][i][kcoor] += f[icoor][jcoor];  // -\delta_{jcoor kcoor} \partial_{icoor} + \delta_{jcoor icoor}\partial_{kcoor}
//                 //                for ( i = 0;i < fe_u.nbNode;i++ )
//                 //                    for(jcoor=0;jcoor<fe_u.nbCoor;++jcoor)
//                 //                        s[i][jcoor] = 0.0;
//                               //for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
//                               //{
//                                 }
//                         }
                        //                    }
                        //                for ( UInt kcoor = 0;kcoor < fe_u.nbCoor;kcoor++ )
                        //                    {
                for ( int kcoor = 0;kcoor < fe_u.nbCoor;kcoor++ )
                    {
                l=0;
                //for ( kcoor = 0;kcoor < fe_u.nbCoor;kcoor++ )
                 for ( icoor = 0;icoor < fe_u.nbCoor;icoor++ )
                 for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
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
    for ( UInt kcoor = 0;(int)kcoor < fe_u.nbCoor;kcoor++ )
        {
            //            for ( short l = 0;i < fe_u.nbNode;i++ )
            //for ( jcoor = 0;jcoor < fe_u.nbCoor;jcoor++ )
            double l = 0.;

            ElemMat::matrix_view mat = elmat.block( iblock, kcoor );
            for ( int j = 0;j < fe_u.nbNode;j++ )
            for ( int i = 0;i < fe_p.nbNode;i++ )
            mat(i,j)=0.;

            // Loop on nodes, i.e. loop on elementary vector components
            for ( int j = 0;j < fe_u.nbNode;j++ )
                for (int i = 0;i < fe_p.nbNode;i++ )
                    {
                        l=0.;
                        // loop on quadrature points
                        for ( ig = 0;ig < fe_u.nbQuadPt;ig++ )
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
//optionally (if full implicit):
//source_mass2
// -rho * ( \grad u^k dw, v  )
//
//convective term
// rho * ( \grad u^k du, v  )
//
void shape_terms_vel( Real rho,
                      Real mu,
                      const ElemVec& uk_loc,
                      const ElemVec& wk_loc,
                      const ElemVec& convect_loc,
                      const ElemVec& pk_loc,
                      ElemMat& elmat,
                      const CurrentFE& fe,
                      const CurrentFE& fe_p,
                      ElemMat& /*elmatP*/,
                      int /*iblock*/,
                      bool fullImplicit,
                      Real alpha,
                      boost::shared_ptr<ElemMat> elmat_convect
                      )
{
    ASSERT_PRE( fe.hasFirstDeriv(),
                "source_mass needs at least the first derivatives" );


    //Real BGrConv[ fe.nbCoor ][ fe.nbCoor ]/*[fe.nbNode]*/;                 // \grad (convect) at a quadrature point
    Real eta[ fe.nbCoor ][ fe.nbCoor ][fe.nbNode][ fe.nbCoor ];                 // I\div d - (\grad d)^T at a quadrature point
    Real uk[ fe.nbQuadPt ][ fe.nbCoor ];              // u^k quadrature points
    Real guk[ fe.nbQuadPt ][ fe.nbCoor ][ fe.nbCoor ];  // \grad u^k at quadrature points
    Real convect[ fe.nbCoor ];                      // convect at quadrature points
    Real convect_eta[ fe.nbQuadPt ][ fe.nbCoor ][fe.nbNode][ fe.nbCoor ];       // (convect)^T [I\div d - (\grad d)^T] at quadrature points
    Real sigma[ fe.nbCoor ][ fe.nbCoor ];             // [-p^k I + 2*mu e(u^k)] a quadrature point
    Real B[ fe.nbQuadPt ][ fe.nbCoor ][ fe.nbCoor ][fe.nbNode][ fe.nbCoor ];  // [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] at each quadrature point
    Real gd[ fe.nbCoor ][ fe.nbCoor ][fe.nbNode][ fe.nbCoor ];                // \grad d at a quadrature point
    Real A2[ fe.nbQuadPt ][ fe.nbCoor ][ fe.nbCoor ][fe.nbNode][ fe.nbCoor ];  // \grad u^k \grad d + [\grad d]^T[\grad u^k]^T  at each quadrature point
    Real aux[ fe.nbQuadPt ][ fe.nbCoor ][fe.nbNode][fe.nbCoor];
    Real BMass[ fe.nbCoor ][ fe.nbCoor ];
    Real auxMass[ fe.nbQuadPt ][fe.nbNode][fe.nbCoor];

    Real s, sB, sA, sG, pk;

    int icoor, jcoor, ig, kcoor, i, j;
    // loop on quadrature points

    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {
        ////INIT
        for(int p=0; p<fe.nbCoor; ++p)
        {
            uk[ig][p]=0.;
            convect[p]=0.;
            for(int q=0; q<fe.nbCoor; ++q)
            {
                guk[ig][p][q]=0.;
                sigma[p][q]=0.;
                BMass[p][q]=0.;
                for(int d=0; d<fe.nbNode; ++d)
                {
                    auxMass[ ig ][d][q]=0.;
                    aux[ig][p][d][q]=0.;
                    convect_eta[ig][p][d][q]=0.;
                    for(int e=0; e<fe.nbCoor; ++e)
                    {
                        gd[p][q][d][e]=0.;
                        eta[p][q][d][e]=0.;
                        A2[ig][p][q][d][e]=0.;
                        B[ig][p][q][d][e]=0.;
                    }
                }
            }
        }
        ////////END INIT

        // loop on space coordindates
        for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
        {

            // each compontent of uk at each quadrature points
            s = 0.0;
            for ( i = 0;i < fe.nbNode;i++ )
                s += fe.phi( (int)i, (int)ig ) * uk_loc.vec() [ i + icoor * fe.nbNode ];
            uk[ ig ][ icoor ] = s;//uk_x(pt_ig), uk_y(pt_ig), uk_z(pt_ig)

            // each compontent of convect at this quadrature point
            s = 0.0;
            for ( i = 0;i < fe.nbNode;i++ )
                s += fe.phi( (int)i, (int)ig ) * convect_loc.vec() [ i + icoor * fe.nbNode ];
            convect[ icoor ] = s;


            // loop  on space coordindates
            for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            {
                sB = 0.0;
                sG = 0.0;
                for ( i = 0;i < fe.nbNode;i++ )
                {
                    gd[ icoor ][ jcoor ][i][jcoor] = fe.phiDer( (int)i, (int)jcoor, (int)ig );

                    sG += fe.phiDer( (int)i, (int)jcoor, (int)ig ) * uk_loc.vec() [ i + icoor * fe.nbNode ]; //  \grad u^k at each quadrature point
                    sB -= fe.phiDer( (int)i, (int)jcoor, (int)ig ) * wk_loc.vec() [ i + icoor * fe.nbNode ]; //  \grad (- w^k) at this quadrature point
                    sA = -fe.phiDer( (int)i, (int)icoor, (int)ig ) /** d_loc.vec() [ i + jcoor * fe.nbNode ]*/; //  - (\grad d) ^T at this quadrature point
                    eta[ icoor ][ jcoor ][i][ jcoor ] = sA; // -(\grad d) ^T at this quadrature point
                }
                guk[ ig ][ icoor ][ jcoor ] = sG; // \grad u^k at each quadrature point
            }
        }

        //!a part of source_stress
        pk = 0.0;
        for ( i = 0;i < fe_p.nbNode;i++ )
            pk += fe_p.phi( (int)i, (int)ig ) * pk_loc.vec() [ i ]; // p^k at this quadrature point

        // sigma = [-p^k I + 2*mu e(u^k)] a quadrature point
        for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
        {
            for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
                sigma[ icoor ][ jcoor ] = mu * ( guk[ig][ icoor ][ jcoor ] + guk[ig][ jcoor ][ icoor ] );
            sigma[ icoor ][ icoor ] -= pk;
        }

        //!building the tensor \f$\eta = [I\nabla\cdot d - (\nabla d)^T]\f$
        Real z[fe.nbCoor];
        for ( i = 0;i < fe.nbNode;i++ )
        {
            //                for ( int kcoor = 0;kcoor < fe.nbCoor;kcoor++ )
            for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            {
                //if(icoor==jcoor)
                z[ jcoor ] = eta[ jcoor ][ jcoor ][i][jcoor];  //! -\delta_{jcoor, kcoor} \partial_{icoor} + \delta_{jcoor ,icoor}\partial_{kcoor}
            }

            for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
            {
                for (kcoor = 0;(int)kcoor < fe.nbCoor;kcoor++ )
                {
                    //if(icoor==jcoor)
                    eta[ jcoor ][ jcoor ][i][kcoor] -= z[ kcoor ];  //! -\delta_{jcoor, kcoor} \partial_{icoor} + \delta_{jcoor, icoor}\partial_{kcoor}
                }
            }

            //!source_mass1

            s = 0;

            for ( kcoor = 0;kcoor < fe.nbCoor;kcoor++ )
            {
                for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                    for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
                        s += BMass[ icoor ][ jcoor ] * eta[ icoor ][ jcoor ][i][kcoor]; // \grad (-w^k):[I\div d - (\grad d)^T] at each quadrature point
                auxMass[ ig ][i][kcoor] = s;
            }


            for ( kcoor = 0;(int)kcoor < fe.nbCoor;kcoor++ )
            {
                for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
                {
                    s = 0.;
                    for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                    {
                        s += convect[ icoor ] * eta[ icoor ][ jcoor ][i][kcoor]; // convect^T [I\div d - (\grad d)^T]
                    }
                    convect_eta[ ig ][ jcoor ][i][kcoor] = s;
                }
            }

            //! At this point we have:
            //!    \f$ v  \nabla u^k \f at each quadrature point
            //!    \f$ v  [I\nabla\cdot d - (\nabla d)^T] \f$ at each quadrature point: convect_A
            //!    \f$ v  \nabla (-w^k):[I\nabla\cdot d - (\nabla d)^T]\f$



            //! source_stress

            for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
            {
                for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
                {
                    for ( kcoor = 0;kcoor < fe.nbCoor;kcoor++ )
                    {
                        s = 0.;
                        for (int zcoor = 0;zcoor < fe.nbCoor;zcoor++ )
                            s += sigma[ icoor ][ zcoor ] * eta[ zcoor ][ jcoor ][i][kcoor];
                        B[ ig ][ icoor ][ jcoor ][i][kcoor] = s;
                    }
                }
            }
            //! source_stress2
            for ( kcoor = 0;kcoor < fe.nbCoor;kcoor++ )
            {
                for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                {
                    for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
                    {
                        s = 0.;
                        for ( int zcoor = 0;zcoor < fe.nbCoor;zcoor++ )
                            // \grad u^k \grad d + [\grad d]^T[\grad u^k]^T  at each quadrature point
                            s += guk[ig][ icoor ][ zcoor ] * gd[ zcoor ][ jcoor ][i][kcoor] + gd[ zcoor ][ icoor ][i][kcoor] * guk[ig][ jcoor ][ zcoor ];
                        A2[ ig ][ icoor ][ jcoor ][i][kcoor] = s;
                    }
                }
            }

            if(fullImplicit)//source_mass2 and convective term derivative
                for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
                {
                    //s = 0.0;
                    for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
                    {
                        aux[ ig ][ icoor ][i][jcoor] = guk[ig][ icoor ][ jcoor ]  * fe.phi( i, ig ) /**alpha*/;
                    }
                }

        }

    }



    //
    // Numerical integration
    //
    Real g=0.;
    // loop on coordinates, i.e. loop on elementary vector blocks
    for ( icoor = 0;icoor < fe.nbCoor;icoor++ )
    {
        for ( kcoor = 0;kcoor < fe.nbCoor;kcoor++ )
        {

            // the block iccor of the elementary vector
            ElemMat::matrix_view mat = elmat.block( icoor, kcoor );
            boost::shared_ptr<ElemMat::matrix_view> mat_convect;

            if(fullImplicit)
                mat_convect.reset(new ElemMat::matrix_view(elmat_convect->block( icoor, kcoor )));
            // loop on nodes, i.e. loop on components of this block
            for ( i = 0;i < fe.nbNode;i++ )
            {
                for ( j = 0;j < fe.nbNode;j++ )
                {

                    // loop on quadrature points
                    s = 0.;
                    g = 0.;

                    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
                    {
                        // - convect^T [I\div d - (\grad d)^T] (u^k)^T :\grad\phi_i
                        for ( jcoor = 0;jcoor < fe.nbCoor;jcoor++ )
                        {
                            //source_mass1
                            s += convect_eta[ ig ][ jcoor ][j][kcoor] * guk[ ig ][ icoor ][jcoor] * fe.phi( (int)i, (int)ig ) * fe.weightDet( ig )*rho;
                            //source_stress
                            s += B[ ig ][ icoor ][ jcoor ][j][kcoor] * fe.phiDer( (int)i, (int)jcoor, (int)ig ) * fe.weightDet( ig );
                            //source_stress2
                            s -= fe.phiDer( (int)i, (int)jcoor, (int)ig ) * A2[ ig ][ icoor ][ jcoor ][j][kcoor] * fe.weightDet( ig )*mu;
                        }
                        if(fullImplicit)
                        {
                            //source_mass1
                            s += auxMass[ ig ][j][kcoor] * uk[ ig ][ icoor ] * fe.phi( i, ig ) * fe.weightDet( ig )*rho;
                            //source_mass2
                            s -= aux[ ig ][ icoor ][j][kcoor] * fe.phi( i, ig ) * fe.weightDet( ig )*rho*alpha;
                            //convective term
                            g += aux[ ig ][ icoor ][j][kcoor] * fe.phi( i, ig ) * fe.weightDet( ig )*rho;
                        }
                    }
                    mat( i , j ) += s;
                    if(fullImplicit)
                        (*mat_convect)( i , j ) += g;
                }
            }
        }
    }



}

}
