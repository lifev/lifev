/* -*- mode: c++ -*-
   This program is part of the LifeV library 
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
#include "darcySolver.hpp"
#include "bc_manage.hpp"
#include "medit_wrtrs.hpp"
#include "clapack.h"
#include "user_diffusion.hpp"

//====================================================

DarcySolver::DarcySolver(const GetPot& data_file):
  DarcyHandler(data_file),
  msrPattern(tpdof),
  mat(msrPattern),
  globalTP(dimTPdof),
  globalF(dimTPdof),
  globalP(dimPdof),
  globalFlux(dimTPdof),
  elvecHyb(refTPFE.nbDof,1),
  elvecSource(pfe.nbNode,1),
  elvecFlux(vfe.nbNode,1),
  elmatMix(vfe.nbNode,1,1, pfe.nbNode,0,1, refTPFE.nbDof,0,1),
  elmatHyb(refTPFE.nbDof,1,1),
  BtB(pfe.nbNode,pfe.nbNode),
  CtC(refTPFE.nbDof, refTPFE.nbDof),
  BtC(pfe.nbNode, refTPFE.nbDof),
  signLocalFace(mesh.numVolumes(),numFacesPerVolume),
  diffusion_scalar_ele(1 /*mesh.numVolumes()*/)  //to save memory when possible...
{

  signLocalFace = -1.;

  // initialisation of the space dependent diffusion with diffusion_scalar
  diffusion_scalar_ele = diffusion_scalar;

  /*
    Initialization of the signs of the faces:
    the rule is the following: if 
  */
  RegionMesh3D<LinearHexa>::VolumeType* vol;
  ID iglobface,iglobvol;
  for(ID ivol = 1; ivol<= mesh.numVolumes(); ivol++){
    vol = &(mesh.volumeList(ivol));
    iglobvol = vol->id();
    for(ID ilocface=1;ilocface<=vol->numLocalFaces;ilocface++){
      iglobface = mesh.localFaceId(iglobvol,ilocface);
      if(mesh.faceElement(iglobface,1) == iglobvol){
	signLocalFace(iglobvol-1,ilocface-1) = 1.;
      }
    }
  }
}


void DarcySolver::_element_computation(int ielem)
{
  // update the current element for the Velocity (only).
  vfe.updatePiola(mesh.volumeList(ielem));
  /*
    update the current element for the Pressure (used for the source term).
    The only necessary part is the quadpt.
    The jacobian is not used. (check...)
  */
  pfe.updateJacQuadPt(mesh.volumeList(ielem));
  /*
    modify the (0,0) block (A) of the matrix
    the blocks (0,1) (B) and (0,2) (C) are independent of the element
    and have already been computed
  */
  // only one block is set to zero since the other ones are not recomputed
  elmatMix.block(0,0) = 0.;

  double xg,yg,zg;

  switch(diffusion_type){
  case 0:
    //-------------------------
    // *** scalar diffusion ***
    //-------------------------
    switch(diffusion_function){
    case 0: // constant diffusion given in the data file
      mass_Hdiv(1./diffusion_scalar,elmatMix,vfe,0,0);
      break;
    case 9: // porous medium with a function
      pfe.barycenter(xg,yg,zg); // coordinate of the barycenter of the current element
      diffusion_scalar = permeability_sd009(xg, yg, zg);
      mass_Hdiv(1./diffusion_scalar,elmatMix,vfe,0,0);
      break;
    case 10: // porous medium with a function
      pfe.barycenter(xg,yg,zg); // coordinate of the barycenter of the current element
      diffusion_scalar = permeability_sd010(xg, yg, zg);
      mass_Hdiv(1./diffusion_scalar,elmatMix,vfe,0,0);
      break;
    default:
      cerr << "Unknown function for scalar diffusion. Change physics/diffusion_function in the data file\n";
      exit(1);
    }
    break;
  case 1:
    //-------------------------
    // *** tensor diffusion ***
    //-------------------------
    {
      KNM<double> permlower(3,3),invpermea(3,3);
      switch(diffusion_function){
      case 0: // constant diffusion tensor given in the data file
	permlower = diffusion_tensor;
	permlower *= diffusion_scalar; 
	/* // the diffusion matrix is divided by the viscosity which sould be 1/diffusion_scalar!
           Remark: not very optimal since in this case
	   this constant matrix is inverted on each elements.
	   (can be easily improved if needed)
	*/
	break;
      case 1: // fibrous medium
	pfe.barycenter(xg,yg,zg); // coordinate of the barycenter of the current element
	permlower = fibrous_permea(1./diffusion_scalar,diffusion_tensor, xg); //! added "1./"
	break;
      default:
	cerr << "Unknown function for tensor diffusion. Change physics/diffusion_function in the data file\n";
	exit(1);
      }
      //
      // we compute the inverse of permlower
      //
      // PERM <- L and Lt where L Lt is the Cholesky factorization of PERM
      int NBT[1] = {3}; //  dim of tensor permeabilite
      int INFO[1] = {0};
      dpotrf_("L", NBT, permlower, NBT, INFO);
      ASSERT_PRE(!INFO[0],"Lapack factorization of PERM is not achieved.");
      dpotri_("L", NBT, permlower, NBT, INFO);
      ASSERT_PRE(!INFO[0],"Lapack solution of PERM is not achieved."); 
      permlower(0,1) = permlower(1,0);
      permlower(0,2) = permlower(2,0);
      permlower(1,2) = permlower(2,1);
      invpermea = permlower;
      //
      mass_Hdiv(invpermea, elmatMix, vfe,0,0); //  modify the (0,0) block of the matrix
      break;
    }
  default:
    cerr << "diffusion_type=" << diffusion_type << " ??? \n";
    exit(1);
  }
}


void DarcySolver::computeHybridMatrixAndSourceRHS()
{
  int INFO[1] = {0};
  int NBRHS[1] = {1};//  nb columns of the rhs := 1.
  int NBP[1] = {pfe.nbNode};//  pressure dof
  int NBU[1] = {vfe.nbNode};//  velocity dof
  int NBL[1] = {refTPFE.nbDof};//  trace of pressure dof (lagrange multiplier)
  double ONE_[1]      = {1.0};
  double MINUSONE_[1] = {-1.0};
  double ZERO_[1]     = {0.0};
  //
  elvecHyb.zero();
  elvecSource.zero();
  elvecFlux.zero();
  elmatMix.zero();
  elmatHyb.zero();
  /*
    Update the divergence Matrix
    (independant of the current element thanks to the Piola transform)
  */
  grad_Hdiv(1.,elmatMix,vfe,pfe,0,1);
  /*
    Update the Boundary Matrix
    (independant of the current element thanks to the Piola transform.)
  */
  TP_VdotN_Hdiv(1., elmatMix, refTPFE, 0, 2 );
  //
  if(verbose>3){
    cerr << "elmatHyb : \n";
    elmatHyb.showMe();
    cerr << "elmatMix : \n";
    elmatMix.showMe();
  }
  //
  globalTP=0.0; 
  globalF=0.0;
  globalP = 0.0;
  globalFlux = 0.0;
  mat.zeros(); 
  //
  //
  /*
  Initial Element Matrix : (where t is transpose)
  
    [ A  B  C    [ U       [ Pext       [ F1      <- F1 : Dirichlet
      Bt 0  0      P    =    F       =    F2      <- F2 : Source term
      Ct 0  0 ]    L ]       UExt ]       F3 ]    <- F3 : Neumann & Robin

   The fluxes U_i are not supposed to be continuous accross the faces: this
   yields a block diagonal matrix A (therefore A^{-1} is not full and can be
   easily computed at the element level).
   
   The L_i are the Lagrange multipliers forcing the continuity of the flux
   across the faces. 
   
   We eliminate U and P and we keep L 
   
    (i) Matrix to be assembled and solved to determine L:
 
    [   Ct A^{-1} C -  Ct A^{-1} B ( Bt A^{-1} B )^{-1} Bt A^{-1} C ] * L
    = 
    Ct A^{-1}  *  F1  -  Ct A^{-1}  B ( Bt A^{-1} B )^{-1} Bt A^{-1}  *  F1  
    +  Ct A^{-1}  B ( Bt A^{-1} B )^{-1} F2
    -  F3
    
    (ii) Recover the pression at the element level:
    
    P = - ( Bt A^{-1} B )^{-1} Bt A^{-1} C * L   
        + ( Bt A^{-1} B )^{-1} Bt A^{-1} * F1 
        - ( Bt A^{-1} B )^{-1}  *  F2

    (iii) Recover the flux at the element level:

    U = -  A^{-1} B  * P  -  A^{-1} C  *  L  +  A^{-1}    * F1

  */
  
  Tab2d AA = elmatMix.block(0,0);
  Tab2d BB = elmatMix.block(0,1);
  Tab2d CC = elmatMix.block(0,2);
  
  SourceAnalyticalFct sourceAnalytical;

  for(UInt ivol = 1; ivol<= mesh.numVolumes(); ivol++){
    //----------------------------
    // LOOP ON THE VOLUME ELEMENTS
    //----------------------------
    _element_computation(ivol);
    //
    AA = elmatMix.block(0,0); // AA <- A
    BB = elmatMix.block(0,1); // BB <- B
    CC = elmatMix.block(0,2); // CC <- C

    //....................
    // MATRIX OPERATIONS
    //....................
    //   AA <- L and Lt where L Lt is the Cholesky factorization of A
    dpotrf_("L", NBU, AA , NBU , INFO );
    ASSERT_PRE(!INFO[0],"Lapack factorization of A is not achieved.");
    // Compute BB <-  L^{-1} B (solve triangular system)
    dtrtrs_("L", "N", "N", NBU, NBP, AA, NBU, BB, NBU, INFO);
    ASSERT_PRE(!INFO[0],"Lapack Computation B = L^{-1} B  is not achieved.");
    // Compute CC <-  L^{-1} C (solve triangular system)
    dtrtrs_("L", "N", "N", NBU, NBL, AA, NBU, CC, NBU, INFO);
    ASSERT_PRE(!INFO[0],"Lapack Computation C = L^{-1} C  is not achieved.");
    // Compute BtB <-  Bt L^{-t} L^{-1} B = Bt A^{-1} B
    // (BtB stored only on lower part)
    dsyrk_("L", "T", NBP, NBU, ONE_, BB, NBU, ZERO_, BtB, NBP);
    // Compute CtC <-  Ct L^{-t} L^{-1} C = Ct A^{-1} C
    // (CtC stored only on lower part)
    dsyrk_("L", "T", NBL, NBU, ONE_, CC, NBU, ZERO_, CtC, NBL);
    // Compute BtC <- Bt L^{-t} L^{-1} C = Bt A^{-1} C
    // (BtC fully stored)
    dgemm_("T", "N", NBP, NBL, NBU, ONE_, BB, NBU, CC, NBU, 
           ZERO_, BtC, NBP );
    //BtB <- LB and LBt where LB LBt is the cholesky factorization of Bt A^{-1} B
    dpotrf_("L", NBP, BtB , NBP, INFO );
    ASSERT_PRE(!INFO[0],"Lapack factorization of BtB is not achieved.");
    // Compute BtC = LB^{-1} BtC <-  LB^{-1} Bt A^{-1} C
    dtrtrs_("L", "N", "N", NBP, NBL, BtB, NBP, BtC, NBP, INFO);
    ASSERT_PRE(!INFO[0],"Lapack Computation BtC = LB^{-1} BtC is not achieved.");
    // CtC = CtC - (BtC)t BtC 
    // (result stored only on lower part)
    // CtC <- Ct A^{-1} C - Ct A^{-t} B (Bt A^{-1} B)^{-1}  Bt A^{-1} C
    dsyrk_("L", "T", NBL, NBP, MINUSONE_, BtC, NBP, ONE_, CtC, NBL);
    //...........................
    // END OF MATRIX OPERATIONS
    //...........................
    /* Sum up:
       AA  <-  L and Lt where L Lt is the Cholesky factorization of A
       BB  <-  L^{-1} B
       CC  <-  L^{-1} C
       BtB <-  LB and LBt where LB LBt is the factorization of Bt A^{-1} B
       BtC <- LB^{-1} Bt A^{-1} C
       CtC <- Ct A^{-1} C - Ct A^{-t} B (Bt A^{-1} B)^{-1}  Bt A^{-1} C
    */
    //....................
    // VECTOR OPERATIONS
    //....................
    // initialize the rhs vectors.    
    elvecSource.zero();     //  source rhs (NBP) 
    elvecHyb.zero();        //  hybrid rhs : inserted in the global vector.
    Tab1dView RHSTP = elvecHyb.block(0);

    // The source term is computed with a test function in the Pressure space.
    switch(test_case){
    case 33: 
      source(sourceAnalytical, elvecSource, pfe, 0);     
    default:
      source(sourceFct, elvecSource, pfe, 0);
    }
    // initialize the rhs vector (clean this some day...)
    Tab1dView rhs = elvecSource.block(0); // corresponds to F2
    // Compute rhs = LB^{-1} rhs <- LB^{-1} F2
    dtrtrs_("L", "N", "N", NBP, NBRHS, BtB, NBP, rhs, NBP, INFO);
    ASSERT_PRE(!INFO[0],"Lapack Computation rhs = LB^{-1} rhs is not achieved.");
    // Compute RHSTP = t(BtC) rhs <- Ct A^{-1} Bt (Bt A^{-1} B)^{-1} F2
    // (fully stored)
    dgemm_("T", "N", NBL, NBRHS, NBP, ONE_, BtC, NBP, rhs, NBP, 
           ZERO_, RHSTP, NBL );
    //........................
    // END OF VECTOR OPERATIONS.
    // (Not finished: for the moment F1 = 0. TO DO : b.c. )
    //........................
    elmatHyb.block(0,0) = CtC;  /* update the hybrid element matrix
				   Everything is stored in the Lower part */
    assemb_mat_symm_lower(mat,elmatHyb,refTPFE,tpdof,
                          mesh.volumeList(ivol).id(),0,0);
    assemb_vec(globalF,elvecHyb,refTPFE,tpdof,
               mesh.volumeList(ivol).id(), 0);
    //-----------------------------------
    // END OF LOOP ON THE VOLUME ELEMENTS
    //-----------------------------------
  }
}

void DarcySolver::applyBC()
{
  bc_manage(mat,globalF,mesh,tpdof,bc,feBd,1.,0.0);
}

void DarcySolver::solveDarcy()
{
  aztecSolveLinearSyst(mat,globalTP.giveVec(),globalF.giveVec(),
		       globalTP.size(),msrPattern);
}

void DarcySolver::computePresFlux()
{
  //! reset to 0
  globalP = 0.0;
  globalFlux = 0.0;

  int INFO[1] = {0};
  int NBRHS[1] = {1};// nb columns of the rhs := 1.
  int INC1[1] = {1};// increment := 1.
  int NBP[1] = {pfe.nbNode};//  pressure dof
  int NBU[1] = {vfe.nbNode};//  velocity dof
  int NBL[1] = {refTPFE.nbDof};//  trace of pressure dof (lagrange multiplier)
  double ONE_[1]      = {1.0};
  double MINUSONE_[1] = {-1.0};
  double ZERO_[1]     = {0.0};

  /*==========================================
    POST PROCESSING
    compute the pressure (Qk or Pk / element) 
    and the velocity (RTk / element) => 2 (opposite) velocities / face
    ==========================================*/

  Vector& global_flux  = globalFlux;
  
  // No need for CtC in this part: only difference. (+ last dsyrk)
  Tab2d AA = elmatMix.block(0,0);
  Tab2d BB = elmatMix.block(0,1);
  Tab2d CC = elmatMix.block(0,2);

  SourceAnalyticalFct sourceAnalytical;

  RegionMesh3D<LinearHexa>::VolumeType* vol;
  ID iglobface;
  for(ID ivol = 1; ivol<= mesh.numVolumes(); ivol++){
    vol = &(mesh.volumeList(ivol));
    _element_computation(ivol);
    // initialize the temporary matrices

    AA = elmatMix.block(0,0);
    BB = elmatMix.block(0,1);
    CC = elmatMix.block(0,2);
    //....................
    // MATRIX OPERATIONS.
    //...................
    //   AA <- L and Lt where L Lt is the Cholesky factorization of A
    dpotrf_("L", NBU, AA , NBU , INFO );
    ASSERT_PRE(!INFO[0],"Lapack factorization of A is not achieved.");
    // Compute BB <-  L^{-1} B (solve triangular system)
    dtrtrs_("L", "N", "N", NBU, NBP, AA, NBU, BB, NBU, INFO);
    ASSERT_PRE(!INFO[0],"Lapack Computation B = L^{-1} B  is not achieved.");
    // Compute CC <-  L^{-1} C (solve triangular system)
    dtrtrs_("L", "N", "N", NBU, NBL, AA, NBU, CC, NBU, INFO);
    ASSERT_PRE(!INFO[0],"Lapack Computation C = L^{-1} C  is not achieved.");
    // Compute BtB <-  Bt L^{-t} L^{-1} B = Bt A^{-1} B
    // (BtB stored only on lower part)
    dsyrk_("L", "T", NBP, NBU, ONE_, BB, NBU, ZERO_, BtB, NBP);
     // Compute BtC <- Bt L^{-t} L^{-1} C = Bt A^{-1} C
    // (BtC fully stored)
    dgemm_("T", "N", NBP, NBL, NBU, ONE_, BB, NBU, CC, NBU, 
           ZERO_, BtC, NBP );
    //BtB <- LB and LBt where LB LBt is the cholesky factorization of Bt A^{-1} B
    dpotrf_("L", NBP, BtB , NBP, INFO );
    ASSERT_PRE(!INFO[0],"Lapack factorization of BtB is not achieved.");
    // Compute BtC = LB^{-1} BtC <-  LB^{-1} Bt A^{-1} C
    dtrtrs_("L", "N", "N", NBP, NBL, BtB, NBP, BtC, NBP, INFO);
    ASSERT_PRE(!INFO[0],"Lapack Computation BtC = LB^{-1} BtC is not achieved.");
    //...........................
    // END OF MATRIX OPERATIONS.
    //..........................
    /* Sum up:
       AA  <-  L and Lt where L Lt is the Cholesky factorization of A
       BB  <-  L^{-1} B
       CC  <-  L^{-1} C
       BtB <-  LB and LBt where LB LBt is the factorization of Bt A^{-1} B
       BtC <- LB^{-1} Bt A^{-1} C
    */
    //...................
    // VECTOR OPERATIONS (Computation of Pressure and Velocities)
    //...................

    //________________________________
    // 1/ Computation of the PRESSURE
    //________________________________
    
    // initialize the rhs vectors.    
    elvecSource.zero();     //  source rhs (NBP) 
    elvecHyb.zero();        //  hybrid rhs : extracted from the global vector.
    // The source term is computed with a test function in the Pressure space.
    // Beware: integrate the source term... ?
    switch(test_case){
    case 33:
      source(sourceAnalytical, elvecSource, pfe, 0);      
    default:
      source(sourceFct, elvecSource, pfe, 0);
    }
    // initialize the rhs vector (clean this some day...)
    Tab1dView rhs = elvecSource.block(0); // corresponds to F2
    // Compute rhs = LB^{-1} rhs <- LB^{-1} F2
    dtrtrs_("L", "N", "N", NBP, NBRHS, BtB, NBP, rhs, NBP, INFO);
    ASSERT_PRE(!INFO[0],
	       "Lapack Computation rhs = LB^{-1} rhs is not achieved.");
    // extract the resulting TP for the current fe and put it into elvecHyb. 
    extract_vec(globalTP, elvecHyb, refTPFE, tpdof,ivol, 0);
    Tab1dView RHSTP = elvecHyb.block(0);
    // RHSTP = elvecHyb.block(0)  contains the local TP for the current fe.
    // rhs = BtC * RHSTP + rhs <- LB^{-1} Bt A^{-1} C * L + LB^{-1} F2 
    /* FAUTE DE SIGNE !!!! version originale @@@@@@@@@@@@@@
      dgemm_("N", "N", NBP, NBRHS, NBL, ONE_, BtC, NBP, RHSTP, NBL, 
      ONE_, rhs, NBP, strlen("T"), strlen("N") );*/
    dgemm_("N", "N", NBP, NBRHS, NBL, MINUSONE_, BtC, NBP, RHSTP, NBL, 
           MINUSONE_, rhs, NBP );
    //  rhs = LB^{-T} rhs
    //  rhs <- - (Bt A^{-1} B)^{-1} Bt A^{-1} C * L - (Bt A^{-1} B)^{-1} F2
     // TO DO : the case when F1 and F3 are not zero
    dtrtrs_("L", "T", "N", NBP, NBRHS, BtB, NBP, rhs, NBP, INFO);
    ASSERT_PRE(!INFO[0],
	       "Lapack Computation rhs = LB^{-T} rhs is not achieved.");
    // contains the pressure for the current element.
    /* Put the pressure of the current fe ("elvecSource")
      in the global vector globalP.*/
    assemb_vec( globalP, elvecSource, refPFE, pdof,ivol, 0);    
    //__________________________________
    // 2/ Computation of the VELOCITIES
    //__________________________________
    // initialize the element flux vector.    
    elvecFlux.zero();     //  Flux (NBU)
    // initialize the flux vector (clean this some day...)
    Tab1dView flux = elvecFlux.block(0);
    flux = elvecFlux.block(0);
    // Compute  flux = BB * rhs <- L^{-1} B P
    dgemv_("N", NBU, NBP, ONE_, BB, NBU, rhs, INC1, ZERO_, flux, INC1 );
    // Compute  flux = - CC * RHSTP - flux <-   - L^{-1} C TP - L^{-1} B P  
    dgemv_("N", NBU, NBL, MINUSONE_, CC, NBL, RHSTP, INC1, MINUSONE_,
	   flux, INC1 );
     // Compute flux = L^{-T} flux <-  - A^{-1} Ct TP - A^{-1} Bt P
    dtrtrs_("L", "T", "N", NBU, NBRHS, AA, NBU, flux, NBU, INFO);
    ASSERT_PRE(!INFO[0],
	       "Lapack Computation flux = L^{-T} flux is not achieved.");
    // TO DO: the case when F1 is not zero (done, no!??)

    //---------------------------------------
    //           BEWARE!!!!
    // THIS IS WRONG IN THE RTk FOR k >= 1 !!!!
    //---------------------------------------

    for(ID ilocface=1;ilocface<=vol->numLocalFaces;ilocface++){
      iglobface = mesh.localFaceId(ivol,ilocface);
      if(mesh.faceElement(iglobface,1) == ivol){
	global_flux[iglobface-1] = flux( ilocface-1 );
      }
    }
    //..........................
    // END OF VECTOR OPERATIONS. 
    //..........................

    //---------------------------------------
    // END OF THE LOOP ON THE VOLUME ELEMENTS
    //---------------------------------------
  }
}

double DarcySolver::computeFluxFlag(int flag)
{
  RegionMesh3D<LinearHexa>::FaceType* face;
  double faceflux,Fl;
  EntityFlag marker;
  Vector& globalFlux_vec = globalFlux;
  //
  //
  Fl = 0.;
  for(ID i=1; i<= mesh.numBFaces(); i++){
    face = &(mesh.faceList(i));
    marker = face->marker();
    faceflux = globalFlux_vec[i-1];
    if(marker==flag){
      //      if(faceflux <= 0) cerr << "Warning: flux <=0 on outlet ??? \n";
      Fl += faceflux;
    }
  }
  return Fl;
}
