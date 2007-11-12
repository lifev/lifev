/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
/* --------------------------------------------------------------------------*
/                                                                            /
/      ...                                                                   /
/                                                                            /
/                                                                            /
/ ASSEMBLY CLASS                                                             /
/                                                                            /
/ #Version 0.1 Experimental:   31/01/2001 Alessandro Veneziani               /
/ #Version 0.2 Added features for vector problems: VBR or global CSR formats /
/                              02/01/2002 Alessandro Veneziani       /=> la specializzazione per il caso simmetrico e' ancora da fare
/ #Version 0.3 6/10/2004 Added features for discontinuous finite elements
/ #Purpose: Container for the assembly routines                              /
/                                                                            /
/                                                                            /
/---------------------------------------------------------------------------*/
#ifndef _ASSEMBLY_GENERIC
#define _ASSEMBLY_GENERIC

#warning: you should use assemb.hpp instead of assembGeneric.hpp when using EpetraMatrices

#include <life/lifefem/currentFE.hpp>
#include <life/lifefem/currentFEDG.hpp>
#include <life/lifefem/currentIFDG.hpp>
#include <life/lifefem/currentBFDG.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/localDofPattern.hpp>
#include <life/lifefem/bcHandler.hpp>

namespace LifeV
{

// For the moment, the selection of the components of the differential operators to be computed is carried out
// inside the differential operator. This is simpler, but not the best: indeed an "if" is processed
// for every integral computation. It would be better to set up a table of simple operators
// and then to plug the correct diff oper component in the by-components loop
// in "assemble". It will be done (I have to understand "how") **** AV January 2002

template <typename Oper, typename DOF, typename RegionMesh,
typename UsrSourceFct, typename Matrix, typename Vector>
void
assemble( Oper oper, const RegionMesh& mesh, CurrentFE& fe,
          const DOF& dof, const UsrSourceFct& source_fct, Matrix& A, Vector& b )
{
    UInt i, ic, jc;
    UInt nc = b.size() / dof.numTotalDof();
    ElemMat elmat( fe.nbNode, nc, nc );
    ElemVec elvec( fe.nbNode, nc );
    //
    for ( i = 1; i <= mesh.numVolumes(); ++i )
    {
        fe.updateFirstDeriv( mesh.volumeList( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();
        for ( ic = 0; ic < nc; ++ic )
        {
            for ( jc = 0; jc < nc; ++jc )
            {
                compute_mat( elmat, oper, fe, ic, jc ); // compute local matrix
                // the previous line would become:  compute_mat(elmat,oper.comp(ic,jc),fe); ****
                assemb_mat( A, elmat, fe, dof, ic, jc ); // assemble local matrix into global one
            }
#ifdef STABILIZED
            compute_vec_stab( source_fct, elvec, fe, ( int ) ic ); // compute local vector
#else

            compute_vec( source_fct, elvec, fe, ( int ) ic ); // compute local vector
#endif
            // the previous line would become:  compute_vec(source_fct.comp(ic),elvec,fe); ****
            assemb_vec( b, elvec, fe, dof, ic ); // assebmle local vector into global one
        }
    }
}
// version with source term depending on time t
template <typename Oper, typename DOF, typename RegionMesh,
typename UsrSourceFct, typename Matrix, typename Vector>
void
assemble( Oper oper, const RegionMesh& mesh, CurrentFE& fe, const DOF& dof,
          const UsrSourceFct& source_fct, Matrix& A, Vector& b, Real const t )
{
    UInt i, ic, jc;
    UInt nc = b.size() / dof.numTotalDof(); // it should be F.Size() if F is for instance a Physical Vector (like the unknown....)
    ElemMat elmat( fe.nbNode, nc, nc );
    ElemVec elvec( fe.nbNode, nc );
    //
    for ( i = 1; i <= mesh.numVolumes(); ++i )
    {
        fe.updateFirstDeriv( mesh.volumeList( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();
        for ( ic = 0; ic < nc; ++ic )
        {
            for ( jc = 0; jc < nc; ++jc )
            {
                compute_mat( elmat, oper, fe, ic, jc ); // compute local matrix
                // the previous line would become:  compute_mat(elmat,oper.comp(ic,jc),fe); ****
                assemb_mat( A, elmat, fe, dof, ic, jc ); // assemble local matrix into global one
            }
#ifdef STABILIZED
            compute_vec_stab( source_fct, elvec, fe, t, ic ); // compute local vector
#else

            compute_vec( source_fct, elvec, fe, t, ic ); // compute local vector
#endif
            // the previous line would become:  compute_vec(source_fct.comp(ic),elvec,fe); ****
            assemb_vec( b, elvec, fe, dof, ic ); // assemble local vector into global one
        }
    }
}
//
//////////////////////////////////
//
//! assembling of operator matrix for mixed FE.
//! Alain, 31/01/02.
//! I have added the number of component of the unknown corresponding
//! to FE1 and FE2: respectively nc1 and nc2. Is it possible to know it
//! through the other parameters ?
template <typename Oper, typename DOF1, typename DOF2, typename RegionMesh,
typename UsrSourceFct, typename Matrix, typename Vector>
void
assemble_mixed( Oper oper, const RegionMesh& mesh, CurrentFE& fe1,
                CurrentFE& fe2, const DOF1& dof1, const DOF2& dof2,
                UInt const nc1, UInt const nc2, const UsrSourceFct& source_fct,
                Matrix& A, Vector& b )
{
    UInt i, ic, jc;
    ElemMat elmat( fe1.nbNode, nc1, 0, fe2.nbNode, 0, nc2 );
    ElemVec elvec( fe1.nbNode, nc1 );
    for ( i = 1; i <= mesh.numVolumes(); ++i )
    {
        fe1.updateFirstDeriv( mesh.volumeList( i ) ); // as updateFirstDer
        fe2.updateFirstDeriv( mesh.volumeList( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();

        for ( ic = 0; ic < nc1; ++ic )
        {
            for ( jc = 0; jc < nc2; ++jc )
            {

                compute_mat_mixed( elmat, oper, fe1, fe2, ic, jc ); // compute local matrix

                assemb_mat_mixed( A, elmat, fe1, fe2, dof1, dof2, ic, jc ); // assemble local matrix into global one
            }
            compute_vec( source_fct, elvec, fe1, ic ); // compute local vector

            assemb_vec( b, elvec, fe1, dof1, ic ); // assemble local vector into global one
        }
    }
}
// version with source term depending on time
template <typename Oper, typename DOF1, typename DOF2, typename RegionMesh,
typename UsrSourceFct, typename Matrix, typename Vector>
void
assemble_mixed( Oper oper, const RegionMesh& mesh, CurrentFE& fe1,
                CurrentFE& fe2, const DOF1& dof1, const DOF2& dof2,
                UInt const nc1, UInt const nc2, const UsrSourceFct& source_fct,
                Matrix& A, Vector& b, Real const t )
{
    UInt i, ic, jc;
    ElemMat elmat( fe1.nbNode, nc1, 0, fe2.nbNode, 0, nc2 );
    ElemVec elvec( fe1.nbNode, nc1 );
    for ( i = 1; i <= mesh.numVolumes(); ++i )
    {
        fe1.updateFirstDeriv( mesh.volumeList( i ) ); // as updateFirstDer
        fe2.updateFirstDeriv( mesh.volumeList( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();

        for ( ic = 0; ic < nc1; ++ic )
        {
            for ( jc = 0; jc < nc2; ++jc )
            {

                compute_mat_mixed( elmat, oper, fe1, fe2, ic, jc ); // compute local matrix

                assemb_mat_mixed( A, elmat, fe1, fe2, dof1, dof2, ic, jc ); // assemble local matrix into global one
            }
            compute_vec( source_fct, elvec, fe1, t, ic ); // compute local vector

            assemb_vec( b, elvec, fe1, dof1, ic ); // assemble local vector into global one
        }
    }
}
///////////////////
// version with source term depending on time t
template <typename Oper, typename DOF, typename RegionMesh,
typename UsrSourceFct, typename Matrix, typename Vector>
void
assemble_symm( Oper oper, const RegionMesh& mesh, CurrentFE& fe, const DOF& dof,
               const UsrSourceFct& source_fct, Matrix& A, Vector& b, Real const t )
{
    UInt i, ic, jc;
    UInt nc = b.size() / dof.numTotalDof(); // it should be F.Size() if F is for instance a Physical Vector (like the unknown....)
    ElemMat elmat( fe.nbNode, nc, nc );
    ElemVec elvec( fe.nbNode, nc );
    //
    for ( i = 1; i <= mesh.numVolumes(); ++i )
    {
        fe.updateFirstDeriv( mesh.volumeList( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();
        for ( ic = 0; ic < nc; ++ic )
        {
            for ( jc = 0; jc < nc; ++jc )
            {
                compute_mat_symm( elmat, oper, fe, ic, jc ); // compute local matrix expoiting symmetry of the operator
                // the previous line would become:  compute_mat(elmat,oper.comp(ic,jc),fe); ****
                assemb_mat( A, elmat, fe, dof, ic, jc ); // assemble local matrix into global one
            }
            compute_vec( source_fct, elvec, fe, t, ic ); // compute local vector
            // the previous line would become:  compute_vec(source_fct.comp(ic),elvec,fe); ****
            assemb_vec( b, elvec, fe, dof, ic ); // assebmle local vector into global one
        }
    }
}
//
//////////////////////////////////
///////////////////
// BLOCK DIAGONAL OPERATOR ASSEMBLING
// It assembles the same scalar operator over the (3) diagonal blocks
template <typename Oper, typename DOF, typename RegionMesh,
typename UsrSourceFct, typename Matrix, typename Vector>
void
assemble_symm_block_diagonal( Oper oper, const RegionMesh& mesh, CurrentFE& fe, const DOF& dof,
                              const UsrSourceFct& source_fct, Matrix& A, Vector& b, Real const t )
{
    UInt i, ic = 0;
    UInt nc = b.size() / dof.numTotalDof(); // it should be F.Size() if F is for instance a Physical Vector (like the unknown....)
    ElemMat elmat( fe.nbNode, nc, nc );
    ElemVec elvec( fe.nbNode, nc );
    //
    for ( i = 1; i <= mesh.numVolumes(); ++i )
    {
        fe.updateFirstDeriv( mesh.volumeList( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();
        compute_mat_symm( elmat, oper, fe ); // compute local matrix expoiting symmetry of the operator

        ic = 0;
        assemb_mat( A, elmat, fe, dof, ic, ic ); // assemble local matrix into global one
        compute_vec( source_fct, elvec, fe, t, ic ); // compute local vector
        assemb_vec( b, elvec, fe, dof, ic ); // assebmle local vector into global one

        for ( ic = 1; ic < nc; ++ic )
        {
            elmat.block( ic, ic ) = elmat.block( 0, 0 );
            assemb_mat( A, elmat, fe, dof, ic, ic ); // assemble local matrix into global one

            compute_vec( source_fct, elvec, fe, t, ic ); // compute local vector
            // the previous line would become:  compute_vec(source_fct.comp(ic),elvec,fe); ****
            assemb_vec( b, elvec, fe, dof, ic ); // assebmle local vector into global one
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
//
// DG matrix assembling
//
////////////////////////////////////////////////////////////////////////////////

// Compute right hand side vector
void compute_vec_DG_BF(const BCHandler& BCh, ElemVec& bfvec, const CurrentBFDG& bfDG, int iblock);

template<typename OperDG, typename OperDGIF, typename OperDGBF,
     typename DOF, typename DOFBYFACE,
     typename RegionMesh, typename UsrSourceFct,
     typename Matrix, typename Vector>
  void assemble_DG(OperDG operDG, OperDGIF operDGIF, OperDGBF operDGBF,
           const RegionMesh& mesh,
           const BCHandler& BCh,
           CurrentFEDG& feDG, CurrentIFDG& feIFDG, CurrentBFDG& feBFDG,
           const DOF& dof, const DOFBYFACE& dofbyface,
           const UsrSourceFct& source_fct,
           Matrix& A, Vector& b)
{
  UInt i, ic, jc;
  UInt nc = b.size() / dof.numTotalDof();
  UInt iAd, iOp;

  ElemMat elmat(feDG.nbNode, nc, nc);
  ElemMat bfmat(feDG.nbNode, nc, nc);
  ElemMat ifmat(2 * feDG.nbNode, nc, nc);

  ElemVec elvec(feDG.nbNode, nc);
  ElemVec bfvec(feDG.nbNode, nc);
  ElemVec ifvec(feDG.nbNode, nc, feDG.nbNode, nc);

/*   ElemMat::matrix_view matIF = ifmat.block(0, 0); */
/*   ElemMat::matrix_view matBF = bfmat.block(0, 0); */

  // Assembling volume integrals contributions
  for(i = 1; i <= mesh.numVolumes(); ++i){
    feDG.updateFirstDeriv(mesh.volumeList(i));

    elmat.zero();
    elvec.zero();

    for(ic = 0; ic < nc; ic++){
      for(jc = 0; jc < nc; jc++){
    compute_mat_DG(elmat, operDG, feDG, ic, jc);
    assemb_mat_DG(A, elmat, feDG, dof, i, ic, jc);
      }
      compute_vec_DG(source_fct, elvec, feDG, ic);
      assemb_vec_DG(b, elvec, feDG, dof, ic);
    }
  }

  // Assembling boundary faces contributions (boundary faces are stored first)
  for(i = 1; i <= mesh.numBFaces(); ++i){
    iAd = (UInt)mesh.faceList(i).ad_first();

    feBFDG.updateBCType(mesh.faceList(i), BCh);
    // The update member to call should be chosen according to the penalization method...
    feBFDG.updateMeasNormalQuadPtFirstDerivAd(mesh.faceList(i), mesh.volumeList(iAd));

    bfmat.zero();
    bfvec.zero();

    for(ic = 0; ic < nc; ic++){
      for(jc = 0; jc < nc; jc++){

    // Compute local matrix
    compute_mat_DG_BF(bfmat, operDGBF, feBFDG, ic, jc);

    // Assemble local matrix into global one
    assemb_mat_DG_BF(A, bfmat, feDG, dof, iAd, ic, jc);

      }
      compute_vec_DG_BF(BCh, bfvec, feBFDG, ic);
      assemb_vec_DG_BF(b, bfvec, feBFDG, iAd, dof, ic);

    }
  }

  // Assembling internal faces contributions
  for(i = mesh.numBFaces() + 1; i <= mesh.numFaces(); ++i){
    iAd = (UInt)mesh.faceList(i).ad_first();
    iOp = (UInt)mesh.faceList(i).ad_second();

    feIFDG.updateMeasNormalQuadPtFirstDerivAd(mesh.faceList(i), mesh.volumeList(iAd));
    feIFDG.updateMeasNormalQuadPtFirstDerivOp(mesh.faceList(i), mesh.volumeList(iOp));

    ifmat.zero();
    ifvec.zero();

    for(ic = 0; ic < nc; ic++){
      for(jc = 0; jc < nc; jc++){
    // Compute local matrix
    compute_mat_DG_IF(ifmat, operDGIF, feDG.nbNode, feIFDG, ic, jc);

    // Assemble local matrix into global one
    assemb_mat_DG_IF(A, ifmat, feDG, dofbyface, i - mesh.numBFaces(), ic, jc);

      }
    }
  }

}

template<typename OperDG,  typename OperDGIF, typename OperDGBF,
     typename DOF, typename DOFBYFACE,
     typename RegionMesh, typename UsrSourceFct,
     typename Matrix, typename Vector, typename Velocity>
  void assemble_AdvecDG(OperDG operDG, OperDGIF operDGIF, OperDGBF operDGBF,
            const RegionMesh& mesh,
            const BCHandler& BCh,
            Velocity& u,
            CurrentFEDG& feDG, CurrentIFDG& feIFDG, CurrentBFDG& feBFDG,
            const DOF& dof, const DOFBYFACE& dofbyface,
            const UsrSourceFct& source_fct,
            Matrix& A, Matrix& M, Vector& b)
{
  UInt i, ic, jc;
  UInt nc = b.size() / dof.numTotalDof();
  UInt iAd, iOp;

  ElemMat elmat(feDG.nbNode, nc, nc);
  ElemMat bfmat(feDG.nbNode, nc, nc);
  ElemMat ifmat(2 * feDG.nbNode, nc, nc);

  ElemVec elvec(feDG.nbNode, nc);
  ElemVec bfvec(feDG.nbNode, nc);
  ElemVec ifvec(feDG.nbNode, nc, feDG.nbNode, nc);

  // Assembling volume integrals contributions
  for(i = 1; i <= mesh.numVolumes(); ++i){
    feDG.updateFirstDerivQuadPtMass(mesh.volumeList(i));

    elmat.zero();
    elvec.zero();

    for(ic = 0; ic < nc; ic++){
      for(jc = 0; jc < nc; jc++){
    compute_mat_DG(elmat, operDG, feDG, ic, jc);
    assemb_mat_DG(A, elmat, feDG, dof, i, ic, jc);

    assemb_mass_DG(M, feDG.mass, feDG, dof, i, ic, jc);
      }
      compute_vec_DG(source_fct, elvec, feDG, ic);
      assemb_vec_DG(b, elvec, feDG, dof, ic);
    }
  }

  // Assembling boundary faces contributions (boundary faces are stored first)
  for(i = 1; i <= mesh.numBFaces(); ++i){
    iAd = (UInt)mesh.faceList(i).ad_first();

    feBFDG.updateBCType(mesh.faceList(i), BCh);
    feBFDG.updateMeasNormalQuadPtFirstDerivAd(mesh.faceList(i), mesh.volumeList(iAd));

    bfmat.zero();
    bfvec.zero();

    for(ic = 0; ic < nc; ic++){
      for(jc = 0; jc < nc; jc++){

    // Compute local matrix
    compute_mat_DG_BF(bfmat, operDGBF, feBFDG, ic, jc);

    // Assemble local matrix into global one
    assemb_mat_DG_BF(A, bfmat, feDG, dof, iAd, ic, jc);
      }
      compute_vec_AdvecDG_BF(BCh, u, bfvec, feBFDG, ic);
      assemb_vec_DG_BF(b, bfvec, feBFDG, iAd, dof, ic);

    }
  }

  // Assembling internal faces contributions
  for(i = mesh.numBFaces() + 1; i <= mesh.numFaces(); ++i){
    iAd = (UInt)mesh.faceList(i).ad_first();
    iOp = (UInt)mesh.faceList(i).ad_second();

    feIFDG.updateMeasNormalQuadPtFirstDerivAd(mesh.faceList(i), mesh.volumeList(iAd));
    feIFDG.updateMeasNormalQuadPtFirstDerivOp(mesh.faceList(i), mesh.volumeList(iOp));

    ifmat.zero();
    ifvec.zero();

    for(ic = 0; ic < nc; ic++){
      for(jc = 0; jc < nc; jc++){
    // Compute local matrix
    compute_mat_DG_IF(ifmat, operDGIF, feDG.nbNode, feIFDG, ic, jc);

    // Assemble local matrix into global one
    assemb_mat_DG_IF(A, ifmat, feDG, dofbyface, i - mesh.numBFaces(), ic, jc);

      }
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////
//
template <typename Oper>
void compute_mat( ElemMat& elmat, Oper& oper,
                  const CurrentFE& fe, int iblock = 0, int jblock = 0 )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt ig;
    unsigned int i, j;
    double s;
    Real x = 0, y = 0, z = 0;
    for ( i = 0;i < ( UInt ) fe.nbNode;i++ )
    {
        for ( j = 0;j < ( UInt ) fe.nbNode;j++ )
        {
            s = 0;
            for ( ig = 0;ig < ( UInt ) fe.nbQuadPt;ig++ )
            {
#ifdef SPACE_DEP_OPERATOR
                fe.coorQuadPt( x, y, z, ig );
#endif

                s += oper( i, j, ig, x, y, z, iblock, jblock ) * fe.weightDet( ig );
                // the previous line would become:  oper(i,j,ig,x,y,z)*fe.weightDet(ig); ****
            }
            mat( ( int ) i, ( int ) j ) += s;
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
//
// Compute local matrices for DG
//
////////////////////////////////////////////////////////////////////////////////


// Volume integrals

template<typename Oper>
void compute_mat_DG(ElemMat& elmat, Oper& oper,
         const CurrentFEDG& fe, int iblock=0,int jblock=0)
{
    ElemMat::matrix_view mat = elmat.block(iblock,jblock);
    UInt ig;
    unsigned int i,j;
    double s;

    Real x=0,y=0,z=0;
    for(i=0;i<(UInt)fe.nbNode;i++){
        for(j=0;j<(UInt)fe.nbNode;j++){
            s = 0.;

            for(ig=0;ig<(UInt)fe.nbQuadPt;ig++){
                fe.coorQuadPt(x,y,z,ig);

                s += oper(i, j, ig, x, y, z, 0, 0, iblock, jblock) * fe.weightDet(ig);
            }

            mat((int)i,(int)j) += s;
        }
    }
}

// Boundary integrals (boundary faces)

template<typename Oper>
void compute_mat_DG_BF(ElemMat& bfmat, Oper& oper,
         const CurrentBFDG& fe, int iblock=0,int jblock=0)
{
  ElemMat::matrix_view mat = bfmat.block(iblock,jblock);
  UInt ig;
  unsigned int i, j;
  double s;
  Real x=0,y=0,z=0;

  for(i=0;i<(UInt)fe.nbNodeAd;i++){
    for(j=0;j<(UInt)fe.nbNodeAd;j++){
      s = 0;
      for(ig=0;ig<(UInt)fe.nbQuadPt;ig++){
    fe.coorQuadPt(x,y,z,ig);

    s += oper(i, j, ig, x, y, z, 0, 0, iblock, jblock) * fe.weightMeas(ig);

       }
      mat((int)i, (int)j) += s;
    }
  }
}

// Boundary integrals (internal faces)

template<typename Oper>
void compute_mat_DG_IF(ElemMat& ifmat, Oper& oper, int nbNodeAd,
         const CurrentIFDG& fe, int iblock=0,int jblock=0)
{
  ElemMat::matrix_view mat = ifmat.block(iblock,jblock);
  UInt ig;
  unsigned int i, j, H, K;
  Real s;

  Real x = 0, y = 0, z = 0;

  for(K = 0; K < 2; K++){
    for(H = 0; H < 2; H++){

      for(i = 0; i < (UInt)fe.nbNodeAd; i++){
    for(j = 0; j < (UInt)fe.nbNodeAd; j++){
      s = 0.;

      for(ig = 0; ig < (UInt)fe.nbQuadPt; ig++){
        fe.coorQuadPt(x,y,z,ig);
        s += oper(i, j, ig, x, y, z, K, H, iblock, jblock) * fe.weightMeas(ig);
      }// for ig

      mat((int)(H * nbNodeAd + i), (int)(K * nbNodeAd + j)) += s;
    }//for j
      }// for i
    }// for H
  }// for K
}

////////////////////////////////////////////////////////////////////////////////
//
//! compute local matrix element for mixed FE
//Alain 1/02/02.
template <typename Oper>
void compute_mat_mixed( ElemMat& elmat, Oper& oper,
                        const CurrentFE& fe1, const CurrentFE& fe2, int iblock = 0, int jblock = 0 )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int ig;
    int i, j;
    double s;
    Real x, y, z;

    for ( i = 0;i < fe1.nbNode;i++ )
    {
        for ( j = 0;j < fe2.nbNode;j++ )
        {
            s = 0;
            // Alain, 08/02/02:
            // REMARK : the choice of quadrature points depends on the definition
            //          of the type FE1 and FE2 in mainAZ.h.
            //          Therefore the same quadrature rule must be chosen
            //          for both FE1 and FE2.
            for ( ig = 0;ig < fe1.nbQuadPt;ig++ )
            {
                fe1.coorQuadPt( x, y, z, ig );
                s += oper( i, j, ig, x, y, z, iblock, jblock ) * fe1.weightDet( ig );
            }
            mat( ( int ) i, ( int ) j ) += s;
        }
    }
}
//
template <typename Oper>
void compute_mat_symm( ElemMat& elmat, Oper& oper,
                       const CurrentFE& fe, int iblock = 0, int jblock = 0 )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt ig;
    unsigned int i, iloc, jloc;
    double s;
    Real x = 0, y = 0, z = 0;
    //
    // diagonal
    //
    for ( i = 0;i < ( UInt ) fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < ( UInt ) fe.nbQuadPt;ig++ )
        {
#ifdef SPACE_DEP_OPERATOR
            fe.coorQuadPt( x, y, z, ig );
#endif

            s += oper( iloc, iloc, ig, x, y, z, iblock, jblock ) * fe.weightDet( ig );
        }
        mat( ( int ) iloc, ( int ) iloc ) += s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < ( UInt ) fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < ( UInt ) fe.nbQuadPt;ig++ )
        {
#ifdef SPACE_DEP_OPERATOR
            fe.coorQuadPt( x, y, z, ig );
#endif

            s += oper( iloc, jloc, ig, x, y, z, iblock, jblock ) * fe.weightDet( ig );
        }

        mat( ( int ) iloc, ( int ) jloc ) += s;
        mat( ( int ) jloc, ( int ) iloc ) += s;
    }
}

//
/////////////////////////////////
//




template <typename DOF, typename Vector>
void
assembleVector( Vector&          vec,
                ElemVec&         elvec,
                const CurrentFE& fe,
                const DOF&       dof,
                int              iblock,
                int              ipos = 0)

{
//    elmat.showMe();
    ElemVec::vector_view vecView = elvec.block( iblock );

    UInt totdof = dof.numTotalDof();

    int i;
    UInt ig;

    UInt eleID = fe.currentLocalId();

    bool verbose = false;

    for ( i = 0 ; i < fe.nbNode ; i++ )
    {
            ig = dof.localToGlobal( eleID, i + 1 ) - 1 + ipos ; //iblock*totdof1;  // damned 1-base vs 0-base !
            vec[ ig + 1 ] += vecView( i );
    }

    if (verbose)
        std::cout << "ok." << std::endl;
}





template <typename DOF, typename Matrix>
void
assembleMatrix( Matrix&          M,
                ElemMat&         elmat,
                const CurrentFE& fe,
                const DOF&       dof,
                int              iblock,
                int              jblock,
                int              ipos,
                int              jpos)

{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    int i, j, k;
    UInt ig, jg;

    UInt eleID = fe.currentLocalId();

    for ( k = 0 ; k < fe.nbPattern ; k++ )
    {
        i = fe.patternFirst ( k );
        j = fe.patternSecond( k );

        if (mat(i,j) != 0.)
            {
                ig = dof.localToGlobal( eleID, i + 1 ) - 1 + ipos ; //iblock*totdof1;  // damned 1-base vs 0-base !
                jg = dof.localToGlobal( eleID, j + 1 ) - 1 + jpos ; //jblock*totdof2;  // damned 1-base vs 0-base !
                M.set_mat_inc( ig, jg, mat( i, j ) );
            }
    }

}





template <typename DOF1, typename DOF2, typename Matrix>
void
assembleMatrix( Matrix&          M,
                ElemMat&         elmat,
                const CurrentFE& fe1,
                const CurrentFE& fe2,
                const DOF1&      dof1,
                const DOF2&      dof2,
                int              iblock,
                int              jblock,
                int              ipos,
                int              jpos,
                bool             verbose = false)

{
//    elmat.showMe();
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    int i, j, k1, k2;
    UInt ig, jg;

    UInt eleID1 = fe1.currentLocalId();
    UInt eleID2 = fe2.currentLocalId();

    for ( k1 = 0 ; k1 < fe1.nbNode ; k1++ )
    {
        i = fe1.patternFirst( k1 );
        for ( k2 = 0 ; k2 < fe2.nbNode ; k2++ )
        {
            j  = fe2.patternSecond( k2 );
            if (mat(i,j) != 0.)
                {
                    ig = dof1.localToGlobal( eleID1, i + 1 ) - 1 + ipos ; //iblock*totdof1;  // damned 1-base vs 0-base !
                    jg = dof2.localToGlobal( eleID2, j + 1 ) - 1 + jpos ; //jblock*totdof2;  // damned 1-base vs 0-base !
                    M.set_mat_inc( ig, jg, mat( i, j ) );
                }
        }
    }

}


template <typename DOF1, typename DOF2, typename Matrix>
void
assembleTransposeMatrix( Matrix&          M,
                         Real             val,
                         ElemMat&         elmat,
                         const CurrentFE& fe1,
                         const CurrentFE& fe2,
                         const DOF1&      dof1,
                         const DOF2&      dof2,
                         int              iblock,
                         int              jblock,
                         int              ipos ,
                         int              jpos )

{
    ElemMat::matrix_view mat = elmat.block( jblock, iblock );

    int i, j, k1, k2;
    UInt ig, jg;

    UInt eleID1 = fe1.currentLocalId();
    UInt eleID2 = fe2.currentLocalId();

    for ( k1 = 0 ; k1 < fe1.nbNode ; k1++ )
    {
        i = fe1.patternFirst( k1 );
        for ( k2 = 0 ; k2 < fe2.nbNode ; k2++ )
        {
            j = fe2.patternSecond( k2 );
            ig = dof1.localToGlobal( eleID1, i + 1 ) - 1 + ipos;  // damned 1-base vs 0-base !
            jg = dof2.localToGlobal( eleID2, j + 1 ) - 1 + jpos;  // damned 1-base vs 0-base !
            M.set_mat_inc( ig, jg, val*mat( j, i ) );
        }
    }
}



//
/////////////////////////////////
//
template <typename Matrix, typename DOF>
void
assemb_mat( Matrix& M, ElemMat& elmat, const CurrentFE& fe, const DOF& dof, int iblock = 0, int jblock = 0 )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt totdof = dof.numTotalDof();
    int i, j, k;
    UInt ig, jg;
    UInt eleId = fe.currentId();
    for ( k = 0 ; k < fe.nbPattern ; k++ )
    {
        i = fe.patternFirst( k );
        j = fe.patternSecond( k );
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
//        std::cout << ig << " " << jg << " " << mat(i, j) << std::endl;

    }
}


//
/////////////////////////////////
//
// Miguel 01/2004: matrix assembling of interior penalty terms
//
template <typename Matrix, typename DOF>
void
assemb_mat( Matrix& M, ElemMat& elmat, const CurrentFE& fe1, const CurrentFE& fe2, const DOF& dof, int iblock = 0, int jblock = 0 )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt totdof = dof.numTotalDof();
    int i, j, k;
    UInt ig, jg;
    UInt eleId1 = fe1.currentId();
    UInt eleId2 = fe2.currentId();
    for ( k = 0 ; k < fe1.nbPattern ; k++ )
    {
        i = fe1.patternFirst( k );
        j = fe2.patternSecond( k );
        ig = dof.localToGlobal( eleId1, i + 1 ) - 1 + iblock * totdof;
        jg = dof.localToGlobal( eleId2, j + 1 ) - 1 + jblock * totdof;
        M.set_mat_inc( ig, jg, mat( i, j ) );
    }
}

//
/////////////////////////////////
//
//! Added by V. Martin 09/2002 (slightly different from the previous function...)
//!    Works with a reference hybrid element, and a given number of the current geo element.
template <typename DOF, typename Matrix>
void assemb_mat( Matrix& M, ElemMat& elmat, const LocalDofPattern& fe, const DOF& dof,
                 const UInt feId, int iblock = 0, int jblock = 0 )
{
    //  if(elmat.nBlockRow()!=1 || elmat.nBlockCol() != 1){
    //    std::cout << "assemble for vector elem mat not yet implemented\n";
    //    exit(1);
    //  }
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt totdof = dof.numTotalDof();
    int i, j, k;
    UInt ig, jg;
    UInt eleId = feId; //! direct use of the fe identity number. (different from the other assemb_mat)

    for ( k = 0 ; k < fe.nbPattern() ; k++ )
    { //! instead of currentFE::nbPattern  ...
        i  = fe.patternFirst( k );
        j  = fe.patternSecond( k );
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
    }
}

////////////////////////////////////////////////////////////////////////////////
//
// Assembling of the elementary matrices for DG
//
////////////////////////////////////////////////////////////////////////////////

template<typename DOF, typename Matrix>
void assemb_mat_DG(Matrix& M, ElemMat& elmat, const CurrentFEDG& fe, const DOF& dof,
        const UInt feId, int iblock = 0,int jblock = 0)
{
  ElemMat::matrix_view mat = elmat.block(iblock, jblock);
  UInt totdof = dof.numTotalDof();
  int i, j, k;
  UInt ig, jg;

  UInt eleId = feId; //! direct use of the fe identity number. (different from the other assemb_mat)

  for(k = 0 ; k < fe.refFE.elPattern.nbPattern() ; k++){

    i = fe.patternFirst(k);
    j = fe.patternSecond(k);

    ig = dof.localToGlobal(eleId, i + 1) - 1 + iblock * totdof;
    jg = dof.localToGlobal(eleId, j + 1) - 1 + jblock * totdof;

    M.set_mat_inc(ig, jg, mat(i, j));
  }
}

template<typename DOF, typename Matrix>
void assemb_mass_DG(Matrix& M, KNM<Real>& mat, const CurrentFEDG& fe, const DOF& dof,
        const UInt feId, int iblock = 0,int jblock = 0)
{
  UInt totdof = dof.numTotalDof();
  int i, j, k;
  UInt ig, jg;

  UInt eleId = feId; //! direct use of the fe identity number. (different from the other assemb_mat)

  for(k = 0 ; k < fe.refFE.elPattern.nbPattern() ; k++){

    i = fe.patternFirst(k);
    j = fe.patternSecond(k);

    ig = dof.localToGlobal(eleId, i + 1) - 1 + iblock * totdof;
    jg = dof.localToGlobal(eleId, j + 1) - 1 + jblock * totdof;

    M.set_mat_inc(ig, jg, mat(i, j));
  }
}

template<typename DOFBYFACE, typename Matrix>
void assemb_mat_DG_IF(Matrix& M, ElemMat& ifmat, const CurrentFEDG& fe, const DOFBYFACE& dofbyface,
        const UInt ifId, int iblock=0,int jblock=0)
{
  ElemMat::matrix_view mat = ifmat.block(iblock, jblock);
  UInt totdof = dofbyface.numTotalDof();

  int i, j, k;
  UInt ig, jg;

  UInt faceId = ifId;

  for(k = 0 ; k < fe.refFE.facePattern.nbPattern(); k++){

    i = fe.patternFirstFaces(k);
    j = fe.patternSecondFaces(k);

    ig = dofbyface.localToGlobal(faceId, i + 1) - 1 + iblock * totdof;
    jg = dofbyface.localToGlobal(faceId, j + 1) - 1 + jblock * totdof;

    M.set_mat_inc(ig, jg, mat(i,j));
  }
}

template<typename DOF, typename Matrix>
void assemb_mat_DG_BF(Matrix& M, ElemMat& bfmat, const CurrentFEDG& fe, const DOF& dof,
        const UInt AdId, int iblock=0,int jblock=0)
{
  ElemMat::matrix_view mat = bfmat.block(iblock, jblock);
  UInt totdof = dof.numTotalDof();

  int i, j, k;
  UInt ig, jg;

  UInt eleId = AdId;

  for(k = 0; k < fe.refFE.elPattern.nbPattern(); k++){

    i = fe.patternFirst(k);
    j = fe.patternSecond(k);

    ig = dof.localToGlobal(eleId, i + 1) - 1 + iblock * totdof;
    jg = dof.localToGlobal(eleId, j + 1) - 1 + jblock * totdof;

    M.set_mat_inc(ig, jg, mat(i, j));
  }
}

////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////
//
//! assembling of the elementary matrix with a symmetric form: use of the LOWER part
//!
//! How to make it more general? For the upper?? (without too much efficiency loss)
//! CAUTION: TESTED ONLY WITH THE STANDARD PATTERN...
//! V. Martin.
template <typename DOF, typename Matrix>
void assemb_mat_symm_lower( Matrix& M, ElemMat& elmat, const LocalDofPattern& fe,
                            const DOF& dof, const UInt feId, int iblock = 0, int jblock = 0 )
{
    //  if(elmat.nBlockRow()!=1 || elmat.nBlockCol() != 1){
    //    std::cout << "assemble for vector elem mat not yet implemented\n";
    //    exit(1);
    //  }
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt totdof = dof.numTotalDof();
    int i, j, k;
    UInt ig, jg;
    UInt eleId = feId; //! direct use of the fe identity number.

    //  std::cout << " Lower case " << std::endl;
    //  std::cout << mat << std::endl;

    //
    // diagonal
    //
    for ( k = 0 ; k < fe.nbDiag() ; k ++ )
    {  //! instead of currentFE::nbDiag
        i = fe.patternFirst( k );
        j = fe.patternSecond( k );
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
        //    std::cout << mat(i,j)  << "\t";
    }
    //
    // extra diagonal : Lower part.
    //
    //  std::cout << "\nextra diag \nnbupper " <<  fe.nbUpper()  << std::endl;
    for ( k = fe.nbDiag() ; k < fe.nbDiag() + fe.nbUpper() ; k ++ )
    {   //! instead of currentFE::nbUpper
        j = fe.patternFirst( k );   //! transpose of upper -> lower CAUTION: TESTED ONLY WITH THE STANDARD PATTERN...
        i = fe.patternSecond( k );  //! transpose of upper -> lower
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
        M.set_mat_inc( jg, ig, mat( i, j ) );
        //    std::cout << mat(i,j)  << "\t";
    }
}

//
/////////////////////////////////
//
//! assembling of the elementary matrix with a symmetric form: use of the UPPER part
//!
//! How to make it more general? For the upper?? (without too much efficiency loss)
//! CAUTION: TESTED ONLY WITH THE STANDARD PATTERN...
//! V. Martin.
template <typename DOF, typename Matrix>
void assemb_mat_symm_upper( Matrix& M, ElemMat& elmat, const LocalDofPattern& fe,
                            const DOF& dof, const UInt feId, int iblock = 0, int jblock = 0 )
{
    //  if(elmat.nBlockRow()!=1 || elmat.nBlockCol() != 1){
    //    std::cout << "assemble for vector elem mat not yet implemented\n";
    //    exit(1);
    //  }
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt totdof = dof.numTotalDof();
    int i, j, k;
    UInt ig, jg;
    UInt eleId = feId; //! direct use of the fe identity number.

    //  std::cout << " Upper case " << std::endl;

    //
    // diagonal
    //
    for ( k = 0 ; k < fe.nbDiag() ; k ++ )
    {  //! instead of currentFE::nbDiag
        i = fe.patternFirst( k );
        j = fe.patternSecond( k );
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
    }
    //
    // extra diagonal : Upper part.
    //
    for ( k = fe.nbDiag() ; k < fe.nbDiag() + fe.nbUpper() ; k ++ )
    {   //! instead of currentFE::nbUpper
        i = fe.patternFirst( k );   //! upper. CAUTION: TESTED ONLY WITH THE STANDARD PATTERN...
        j = fe.patternSecond( k );  //! upper
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
        M.set_mat_inc( jg, ig, mat( i, j ) );
    }
}

//
//////////////////
//
//! assembling of the elementary matrix into the global one for
//! mixed FE
//Alain 1/02/02.
template <typename DOF1, typename DOF2, typename Matrix>
void
assemb_mat_mixed( Matrix& M, ElemMat& elmat, const CurrentFE& fe1, const CurrentFE& fe2,
                  const DOF1& dof1, const DOF2& dof2, int iblock = 0, int jblock = 0 )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt totdof1 = dof1.numTotalDof();
    UInt totdof2 = dof2.numTotalDof();
    int i, j, k1, k2;
    UInt ig, jg;
    UInt eleID1 = fe1.currentId();
    UInt eleID2 = fe2.currentId();

    for ( k1 = 0 ; k1 < fe1.nbNode ; k1++ )
    {
        i = fe1.patternFirst( k1 );
        for ( k2 = 0 ; k2 < fe2.nbNode ; k2++ )
        {
            j = fe2.patternSecond( k2 );
            ig = dof1.localToGlobal( eleID1, i + 1 ) - 1 + iblock * totdof1;  // damned 1-base vs 0-base !
            jg = dof2.localToGlobal( eleID2, j + 1 ) - 1 + jblock * totdof2;  // damned 1-base vs 0-base !
            M.set_mat_inc( ig, jg, mat( i, j ) );
        }
    }
}

/*!
  Assemble the transposed of an element matrix (multiplied by a
  constant "mulfac") JFG 03/09/02.
*/
template <typename DOF1, typename DOF2, typename Matrix>
void
assemb_tr_mat_mixed( Real mulfac, Matrix& M, ElemMat& elmat,
                     const CurrentFE& fe1, const CurrentFE& fe2,
                     const DOF1& dof1, const DOF2& dof2,
                     int iblock = 0, int jblock = 0 )
{
    ElemMat::matrix_view mat = elmat.block( jblock, iblock );
    UInt totdof1 = dof1.numTotalDof();
    UInt totdof2 = dof2.numTotalDof();
    int i, j, k1, k2;
    UInt ig, jg;
    UInt eleID1 = fe1.currentId();
    UInt eleID2 = fe2.currentId();

    for ( k1 = 0 ; k1 < fe1.nbNode ; k1++ )
    {
        i = fe1.patternFirst( k1 );
        for ( k2 = 0 ; k2 < fe2.nbNode ; k2++ )
        {
            j = fe2.patternSecond( k2 );
            ig = dof1.localToGlobal( eleID1, i + 1 ) - 1 + iblock * totdof1;  // damned 1-base vs 0-base !
            jg = dof2.localToGlobal( eleID2, j + 1 ) - 1 + jblock * totdof2;  // damned 1-base vs 0-base !
            M.set_mat_inc( ig, jg, mulfac * mat( j, i ) ); // note the mat(j,i) which is the
            // only difference with assemb_mat_mixed
        }
    }
}
//
//////////////////
//
void compute_vec( Real constant, ElemVec& elvec, const CurrentFE& fe, int iblock = 0 );
//////////////////
//
// Suggestion by Luca F.
/*
  Why dont we let CurrentFE to store the time as well?
  we may choose then a common layout for the user functions that always has
  the time as entry point.
 */
//
//
//////////////////
//
template <typename UsrFct>
void compute_vec( const UsrFct& fct, ElemVec& elvec, const CurrentFE& fe, int iblock = 0 )
{
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s, x, y, z;
    for ( i = 0;i < fe.nbNode;i++ )
    {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            fe.coorQuadPt( x, y, z, ig );
            s += fe.phi( i, ig ) * fct( x, y, z, iblock + 1 ) * fe.weightDet( ig );
            //the previous line would become: s += fe.phi(i,ig) * fct(x,y,z) * fe.weightDet(ig); ****
        }
        vec( i ) += s;
    }
}
//
template <typename UsrFct>
void compute_vec( const UsrFct& fct, ElemVec& elvec, const CurrentFE& fe )
{
    compute_vec( fct, elvec, fe, 0 );
}
////////////////////////////////////////////////////////////////////////////////
//
// rhs computation for DG elements
//
///////////////////////////////////////////////////////////////////////////////

template<typename UsrFct>
void compute_vec_DG(const UsrFct& fct, ElemVec& elvec, const CurrentFEDG& feDG, int iblock){
  int i,ig;
  ElemVec::vector_view vec = elvec.block(iblock);
  Real s, x, y, z;
  for(i = 0; i < feDG.nbNode; i++){
    s = 0.;
    for(ig = 0; ig < feDG.nbQuadPt; ig++){
      feDG.coorQuadPt(x, y, z, ig);
      s += feDG.phi(i,ig) * fct(x,y,z,iblock) * feDG.weightDet(ig);
    }
    vec(i) += s;
  }
}


// Compute right hand side vector for pure hyperbolic problems
template<typename Velocity>
void compute_vec_AdvecDG_BF(const BCHandler& BCh, Velocity& u, ElemVec& bfvec, const CurrentBFDG& bfDG, int iblock){

  int i, ig, icoor;

  ElemVec::vector_view vec = bfvec.block(iblock);

  Real s, x, y, z, u_normal;

  const BCBase& CurrBC = BCh.GetBCWithFlag(bfDG.marker);
  ASSERT_PRE(bfDG.bcType == Essential, "Only essential conditions admitted in pure hyperbolic problems")

  for(i = 0; i < bfDG.nbNodeAd; i++){
    s = 0.;
    u_normal = 0.;

    for(ig = 0; ig < bfDG.nbQuadPt; ig++){
      bfDG.coorQuadPt(x, y, z, ig);

      for(icoor = 0; icoor < bfDG.nbCoorAd; icoor++){
        u_normal += bfDG.normal(icoor, ig) * u(x, y, z, icoor);
      }

      if(u_normal < 0){
    for(icoor = 0; icoor < bfDG.nbCoorAd; icoor++){
      s += - u_normal * bfDG.phiAd(i, ig) * CurrBC(0., x, y, z, iblock) * bfDG.weightMeas(ig);
    }
      }
    }

    vec(i) += s;
  }
}

// ! versione per il caso stabilizzato (o comunque son termine forzante wrappato)
template <typename OperFct>
void compute_vec_stab( OperFct& fct, ElemVec& elvec, const CurrentFE& fe,
                       int iblock )
{
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s, x, y, z;
    for ( i = 0;i < fe.nbNode;i++ )
    {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            fe.coorQuadPt( x, y, z, ig );
            s += fct( i, ig, x, y, z, iblock + 1 ) * fe.weightDet( ig );
        }
        vec( i ) += s;
    }
}

// version with source term depending on time
template <typename UsrFct>
void compute_vec( const UsrFct& fct, ElemVec& elvec, const CurrentFE& fe, const Real& t, int iblock = 0 )
{
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s, x, y, z;
    for ( i = 0;i < fe.nbNode;i++ )
    {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            x = fe.quadPt( ig, 0 );
            y = fe.quadPt( ig, 1 );
            z = fe.quadPt( ig, 2 );
            s += fe.phi( i, ig ) * fct( t, x, y, z, iblock + 1 ) * fe.weightDet( ig );
            //the previous line would become: s += fe.phi(i,ig) * fct(x,y,z,t) * fe.weightDet(ig); ****
        }
        vec( i ) += s;
    }
}
// ! Stabilized case with time dependence
template <typename OperFct>
void compute_vec_stab( OperFct& fct, ElemVec& elvec, const CurrentFE& fe, Real t,
                       int iblock )
{
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s, x, y, z;
    for ( i = 0;i < fe.nbNode;i++ )
    {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            fe.coorQuadPt( x, y, z, ig );
            s += fct( i, ig, x, y, z, t, iblock + 1 ) * fe.weightDet( ig );
            //the previous line would become: s += fe.phi(i,ig) * fct(x,y,z) * fe.weightDet(ig); ****
        }
        vec( i ) += s;
    }
}
//


template <typename DOF, typename Vector, typename ElemVec>
void
assemb_vec( Vector& V, ElemVec& elvec, const CurrentFE& fe, const DOF& dof, int iblock=0 )
{
    UInt totdof = dof.numTotalDof();
    typename ElemVec::vector_view vec = elvec.block( iblock );
    int i;
    UInt ig;
    UInt eleId = fe.currentId();
    for ( i = 0 ; i < fe.nbNode ; i++ )
    {
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;
        V[ ig ] += vec( i );
    }
}
///

////////////////////////////////////////////////////////////////////////////////
//
// Added by D. A. Di Pietro: vector assembly for DG
//
////////////////////////////////////////////////////////////////////////////////

template<typename DOF, typename Vector, typename ElemVec>
void assemb_vec_DG(Vector& V, ElemVec& elvec,const CurrentFEDG& feDG, const DOF& dof, int iblock)
{
  UInt totdof = dof.numTotalDof();
  typename ElemVec::vector_view vec = elvec.block(iblock);
  int i;
  UInt ig;
  UInt eleId = feDG.currentId();
  for(i = 0; i < feDG.nbNode; i++){
    ig = dof.localToGlobal(eleId, i + 1) - 1 + iblock * totdof;

    V[ig] += vec(i);
  }
}
template<typename DOF, typename Vector, typename ElemVec>
void assemb_vec_DG_BF(Vector& V,ElemVec& bfvec, const CurrentBFDG& bfDG, UInt iAd, const DOF& dof, int iblock){
  UInt totdof = dof.numTotalDof();
  typename ElemVec::vector_view vec = bfvec.block(iblock);
  int i;
  UInt ig;

  UInt eleId = iAd;

  for(i = 0; i < bfDG.nbNodeAd; i++){
    ig = dof.localToGlobal(eleId, i + 1) - 1 + iblock * totdof;

    V[ig] += vec(i);
  }
}

////////////////////////////////////////////////////////////////////////////////

///
//////////////////
/// V. Martin  09/2002
//! version of assemb_vec that works with a LocalDofPattern (and also a RefHybridFE)...
template <typename DOF, typename Vector, typename ElemVec>
void
assemb_vec( Vector& V, ElemVec& elvec, const LocalDofPattern& fe, const DOF& dof,
            const UInt feId, int iblock )
{
    //  if(elvec.nBlockRow()!=1){
    //    std::cout << "assemble for vector elem vec not yet implemented\n";
    //    exit(1);
    //  }
    UInt totdof = dof.numTotalDof();
    typename ElemVec::vector_view vec = elvec.block( iblock );
    int i;
    //  std::cout << "in assemb_vec" << std::endl;
    UInt ig;
    //  UInt eleId = feId;  //simplify.
    for ( i = 0 ; i < fe.nbLocalDof ; i++ )
    {    //! instead of CurrentFE::nbNode
        ig = dof.localToGlobal( feId, i + 1 ) - 1 + iblock * totdof;
        //    std::cout << "i= " << i << std::endl;
        //    std::cout << "ig= " << ig << std::endl;
        V[ ig ] += vec( i );
    }
}
///

///
//////////////////
/// V. Martin  09/2002
//! \function extract_vec
//! \brief from a global vector, extract a vector corresponding to the element feId.
//! works with a LocalDofPattern (and also a RefHybridFE)...
template <typename DOF, typename Vector, typename ElemVec>
void
extract_vec( Vector& V, ElemVec& elvec, const LocalDofPattern& fe, const DOF& dof,
             const UInt feId, int iblock )
{
    //  if(elvec.nBlockRow()!=1){
    //    std::cout << "assemble for vector elem vec not yet implemented\n";
    //    exit(1);
    //  }

    UInt totdof = dof.numTotalDof();
    typename ElemVec::vector_view vec = elvec.block( iblock );
    int i;
    //  std::cout << "in assemb_vec" << std::endl;
    UInt ig;
    for ( i = 0 ; i < fe.nbLocalDof ; i++ )
    {
        ig = dof.localToGlobal( feId, i + 1 ) - 1 + iblock * totdof;
        //    std::cout << "i= " << i << std::endl;
        //    std::cout << "ig= " << ig << std::endl;
        vec( i ) = V[ ig ];
    }
}
}

#endif
