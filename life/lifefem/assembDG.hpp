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
    @file Assemb for DG method.
    @brief This file includes the assemb functions for the use with DG

    @author Simone Deparis <simone.deparis@epfl.ch>
    @date 09 Mar 2010

    The contents of assemb regarding DG have been moved here
 */

#ifndef ASSEMBDG_H
#define ASSEMBDG_H 1

#include <life/lifefem/assemb.hpp>

#include <life/lifefem/currentFEDG.hpp>
#include <life/lifefem/currentIFDG.hpp>
#include <life/lifefem/currentBFDG.hpp>

namespace LifeV {

////////////////////////////////////////////////////////////////////////////////
//
// DG matrix assembling
//
////////////////////////////////////////////////////////////////////////////////

// Compute right hand side vector
void compute_vec_DG_BF(const BCHandler& BCh, ElemVec& bfvec, const CurrentBFDG& bfDG, int iblock);

template<typename OperDG, typename OperDGIF, typename OperDGBF,
     typename DOF, typename DOFBYFACE,
     typename RegionMesh, typename UsrSourceFct>
  void assemble_DG(OperDG operDG, OperDGIF operDGIF, OperDGBF operDGBF,
           const RegionMesh& mesh,
           const BCHandler& BCh,
           CurrentFEDG& feDG, CurrentIFDG& feIFDG, CurrentBFDG& feBFDG,
           const DOF& dof, const DOFBYFACE& dofbyface,
           const UsrSourceFct& source_fct,
           EpetraMatrix<double>& A, EpetraVector& b)
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
  for(i = 1; i <= mesh.numElements(); ++i){
    feDG.updateFirstDeriv(mesh.element(i));

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
    feBFDG.updateMeasNormalQuadPtFirstDerivAd(mesh.faceList(i), mesh.element(iAd));

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

    feIFDG.updateMeasNormalQuadPtFirstDerivAd(mesh.faceList(i), mesh.element(iAd));
    feIFDG.updateMeasNormalQuadPtFirstDerivOp(mesh.faceList(i), mesh.element(iOp));

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
         typename Velocity>
  void assemble_AdvecDG(OperDG operDG, OperDGIF operDGIF, OperDGBF operDGBF,
            const RegionMesh& mesh,
            const BCHandler& BCh,
            Velocity& u,
            CurrentFEDG& feDG, CurrentIFDG& feIFDG, CurrentBFDG& feBFDG,
            const DOF& dof, const DOFBYFACE& dofbyface,
            const UsrSourceFct& source_fct,
            EpetraMatrix<double>& A, EpetraMatrix<double>& M, EpetraVector& b)
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
  for(i = 1; i <= mesh.numElements(); ++i){
    feDG.updateFirstDerivQuadPtMass(mesh.element(i));

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
    feBFDG.updateMeasNormalQuadPtFirstDerivAd(mesh.faceList(i), mesh.element(iAd));

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

    feIFDG.updateMeasNormalQuadPtFirstDerivAd(mesh.faceList(i), mesh.element(iAd));
    feIFDG.updateMeasNormalQuadPtFirstDerivOp(mesh.faceList(i), mesh.element(iOp));

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
// Assembling of the elementary matrices for DG
//
////////////////////////////////////////////////////////////////////////////////

template<typename DOF>
void assemb_mat_DG(EpetraMatrix<double>& M, ElemMat& elmat, const CurrentFEDG& fe, const DOF& dof,
        const UInt feId, int iblock = 0,int jblock = 0)
{
  ElemMat::matrix_view mat = elmat.block(iblock, jblock);
  UInt totdof = dof.numTotalDof();
  return assembleMatrix( M, feId, feId, mat, fe, fe, dof, dof, iblock * totdof, jblock * totdof);

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

template<typename DOF>
void assemb_mass_DG(EpetraMatrix<double>& M, KNM<Real>& mat, const CurrentFEDG& fe, const DOF& dof,
        const UInt feId, int iblock = 0,int jblock = 0)
{


  UInt totdof = dof.numTotalDof();
  return assembleMatrix( M, feId, feId, mat, fe, fe, dof, dof, iblock * totdof, jblock * totdof);

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

template<typename DOFBYFACE>
void assemb_mat_DG_IF(EpetraMatrix<double>& M, ElemMat& ifmat, const CurrentFEDG& fe, const DOFBYFACE& dofbyface,
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

template<typename DOF>
void assemb_mat_DG_BF(EpetraMatrix<double>& M, ElemMat& bfmat, const CurrentFEDG& fe, const DOF& dof,
        const UInt AdId, int iblock=0,int jblock=0)
{
  ElemMat::matrix_view mat = bfmat.block(iblock, jblock);
  UInt totdof = dof.numTotalDof();
  return assembleMatrix( M, AdId, AdId, mat, fe, fe, dof, dof, iblock * totdof, jblock * totdof);

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

////////////////////////////////////////////////////////////////////////////////
//
// Added by D. A. Di Pietro: vector assembly for DG
//
////////////////////////////////////////////////////////////////////////////////

template<typename DOF, typename ElemVec>
void assemb_vec_DG(EpetraVector& V, ElemVec& elvec,const CurrentFEDG& feDG, const DOF& dof, int iblock)
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
template<typename DOF, typename ElemVec>
void assemb_vec_DG_BF(EpetraVector& V,ElemVec& bfvec, const CurrentBFDG& bfDG, UInt iAd, const DOF& dof, int iblock){
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



} // Namespace LifeV

#endif /* ASSEMBDG_H */
