/* -*- mode: c++ -*-
   This program is part of the LifeV library

   Author(s): Vincent Martin <vincent.martin@mate.polimi.it>
   Date: 2004-10-11

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
#ifndef _DARCY_SOLVER_H
#define _DARCY_SOLVER_H

#include "darcyHandler.hpp"
#include "bcManage.hpp"

#include "clapack.h"
#include "user_diffusion.hpp"
#include "medit_wrtrs.hpp"
#include "ensight7Writer.hpp"
#include "sobolevNorms.hpp"


namespace LifeV
{
/*!
  \brief A mixed hybrid Darcy solver
  \file darcySolver.hpp
  \author J.-F. Gerbeau and V. Martin
  \date 11/2002

  Templated class to be used both for tetra or hexa meshes.
  can be used as follows:

  //! HEXA
  DarcySolver< RegionMesh3D<LinearHexa> >
  DarcySolverHexaRT0( data_file, feHexaRT0, feHexaQ0, feHexaRT0Hyb,
                      feHexaRT0VdotNHyb, feHexaQ1,
                      quadRuleHexa8pt, quadRuleQuad4pt );
  //! TETRA
  DarcySolver< RegionMesh3D<LinearTetra> >
  DarcySolverTetraRT0( data_file, feTetraRT0, feTetraP0, feTetraRT0Hyb,
                       feTetraRT0VdotNHyb, feTetraP1,
                       quadRuleTetra15pt, quadRuleTria4pt );

*/

template <typename Mesh>
class DarcySolver:
        public DarcyHandler<Mesh>
{
public:
    //! Constructor
    /*!
      \param data_file GetPot data file
      \param refFE_u reference FE for the RTk velocity
      \param refFE_p reference FE for the Pk/Qk pressure
      \param refFE_tp reference FE for the "hybrid" trace of pressure
      \param refFE_vdotn reference FE for the velocity on the faces
      \param refFE_pnodal reference FE for the P1/Q1 pressure (used for post process)
      \param Qr_u volumic quadrature rule
      \param bdQr_u surface quadrature rule
    */
    DarcySolver( const GetPot& data_file, const RefHdivFE& refFE_u,
                 const RefFE& refFE_p, const RefHybridFE& refFE_tp,
                 const RefHybridFE& refFE_vdotn, const RefFE& refFE_pnodal,
                 const QuadRule& qr_u, const QuadRule& bdqr_u );

private:
    MSRPatt msrPattern;
    MSRMatr<double> mat;
    ScalUnknown<Vector> globalTP;  //!< Trace of Pressure (TP)
    ScalUnknown<Vector> globalF;   //!< Source term
    ScalUnknown<Vector> globalP;   //!< Pressure (P)
    ScalUnknown<Vector> globalFlux;//!< Flux (U)

    ElemVec elvecHyb;//!< Element vector for the Rhs (lives in RTkHyb(fe))
    ElemVec elvecSource;//!< Source element vector (lives in Qk(fe))
    ElemVec elvecFlux;  //!< Flux element vector (lives in RTk(fe))
    ElemMat elmatMix;//!< Mixed hybrid Matrix: [A B C] in [A B C; Bt 0 0; Ct 0 0]
    ElemMat elmatHyb;//!< Hybridization Matrix

    KNM<Real> BtB; //!< tmp array (NBP x NBP)
    KNM<Real> CtC; //!< tmp array (NBL x NBL): storage of the final hybrid array.
    KNM<Real> BtC; //!< tmp array (NBP x NBL)
    /*!
      signLocalFace(ivol,ilocface)
      =  1  if ivol is the first adjacent of face ilocface
      = -1                 second

      Remark: this array could be in regionMesh.
    */
    KNM<Real> signLocalFace;
    KN<Real> diffusion_scalar_ele; //! scalar diffusion coeff, element by element

    SourceFct sourceFct;

private:
    void _element_computation(int i); //!< computations of element matrices

public:
    void computeHybridMatrixAndSourceRHS(); //!< compute the matrix for TP
    //!< and the contribution from the source term to the globalF right hand side.

    void applyBC(); //!< apply the b.c. for the TP problem
    void solveDarcy();//!< solve the linear system for TP with Aztec
    void computePresFlux();//!< Compute P and U (once TP known)
    void postProcessTraceOfPressureRT0();//!< postprocess TP constant per face
    void postProcessVelocityRT0();//!< postprocess Velocity (RT0 per element)
    void postProcessPressureQ0();//!< postprocess P constant per element
    void projectPressureQ1( ScalUnknown<Vector>& nodalPres ); //!< projection of P (Q0/P0) on Q1/P1.
    void postProcessPressureQ1(); //!< postproc of Q1/P1 pressure.
    void projectVelocityQ1( PhysVectUnknown<Vector>& nodalVel ); //!< projection of U (RT0) on Q1/P1.
    void postProcessVelocityQ1();//!<  postproc of Q1/P1 velocity.
    void postProcessEnsight(); //!< postprocessing in ensight format of P and U
    Real computeFluxFlag(int flag);
};


// --------------------------------------------------
// IMPLEMENTATION
// --------------------------------------------------

// constructor
template <typename Mesh>
DarcySolver<Mesh>::DarcySolver( const GetPot& data_file, const RefHdivFE& refFE_u,
                                const RefFE& refFE_p, const RefHybridFE& refFE_tp,
                                const RefHybridFE& refFE_vdotn, const RefFE& refFE_pnodal,
                                const QuadRule& qr_u, const QuadRule& bdqr_u):
    DarcyHandler<Mesh>( data_file, refFE_u, refFE_p, refFE_tp,
                        refFE_vdotn, refFE_pnodal, qr_u, bdqr_u ),
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
    signLocalFace(this->_mesh.numVolumes(),numFacesPerVolume),
    diffusion_scalar_ele(1 /*this->_mesh.numVolumes()*/)  //to save memory when possible...
{

    signLocalFace = -1.;

    // initialisation of the space dependent diffusion with diffusion_scalar
    diffusion_scalar_ele = diffusion_scalar;

    /*
      Initialization of the signs of the faces:
      the rule is the following: if
    */
    ID iglobface,iglobvol;
    for( ID ivol = 1 ; ivol <= this->_mesh.numVolumes() ; ivol++ ) {
        iglobvol = this->_mesh.volumeList(ivol).id();
        for( ID ilocface=1 ; 
             ilocface <= this->_mesh.volumeList(ivol).numLocalFaces ; 
             ilocface++ ) {

            iglobface = this->_mesh.localFaceId(iglobvol,ilocface);
            if( this->_mesh.faceElement(iglobface,1) == iglobvol ) {
                signLocalFace(iglobvol-1,ilocface-1) = 1.;
            }
        }
    }
}


template <typename Mesh>
void DarcySolver<Mesh>::_element_computation(int ielem)
{
    // update the current element for the Velocity (only).
    vfe.updatePiola(this->_mesh.volumeList(ielem));
    /*
      update the current element for the Pressure (used for the source term).
      The only necessary part is the quadpt.
      The jacobian is not used. (check...)
    */
    pfe.updateJacQuadPt(this->_mesh.volumeList(ielem));
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
            std::cerr << "Unknown function for scalar diffusion. Change physics/diffusion_function in the data file\n"
                      << std::endl;
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
                std::cerr << "Unknown function for tensor diffusion. Change physics/diffusion_function in the data file\n"
                          << std::endl;
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
        std::cerr << "diffusion_type=" << diffusion_type << " ??? \n" 
                  << std::endl;
        exit(1);
    }
}

template <typename Mesh>
void DarcySolver<Mesh>::computeHybridMatrixAndSourceRHS()
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
    TP_VdotN_Hdiv(1., elmatMix, refTPFE,refVdotNFE, 0, 2 );
    //
    if(verbose>3){
        std::cout << "elmatHyb : \n" << std::endl;
        elmatHyb.showMe();
        std::cout << "elmatMix : \n" << std::endl;
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

    for(UInt ivol = 1; ivol<= this->_mesh.numVolumes(); ivol++){
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
                              this->_mesh.volumeList(ivol).id(),0,0);
        assemb_vec(globalF,elvecHyb,refTPFE,tpdof,
                   this->_mesh.volumeList(ivol).id(), 0);
        //-----------------------------------
        // END OF LOOP ON THE VOLUME ELEMENTS
        //-----------------------------------
    }
}

template <typename Mesh>
void DarcySolver<Mesh>::applyBC()
{
    bcManage(mat,globalF,this->_mesh,tpdof,bc,feBd,1.,0.0);
}

template <typename Mesh>
void DarcySolver<Mesh>::solveDarcy()
{
    aztecSolveLinearSyst(mat,globalTP.giveVec(),globalF.giveVec(),
                         globalTP.size(),msrPattern);
}

template <typename Mesh>
void DarcySolver<Mesh>::computePresFlux()
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

    ID iglobface;
    for( ID ivol = 1 ; ivol <= this->_mesh.numVolumes(); ivol++ ) {
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
        // rhs contains the pressure for the current element.
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

        for( ID ilocface=1 ;
             ilocface <= this->_mesh.volumeList(ivol).numLocalFaces ; 
             ilocface++ ) {
            iglobface = this->_mesh.localFaceId( ivol , ilocface );
            if( this->_mesh.faceElement( iglobface , 1 ) == ivol ){
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


//-------------------------------------------------
//! post processing part.
//-------------------------------------------------
//! switch to compute or not the analytical solution
#define ANALYTICAL_SOL 0

template <typename Mesh>
Real DarcySolver<Mesh>::computeFluxFlag(int flag)
{
    double faceflux,Fl;
    EntityFlag marker;
    Vector& globalFlux_vec = globalFlux;
    //
    //
    Fl = 0.;
    for( ID iface = 1 ; iface <= this->_mesh.numBFaces() ; iface++ ) {
        marker = this->_mesh.faceList(iface).marker();
        faceflux = globalFlux_vec[ iface - 1 ];
        if( marker == flag ) {
            Fl += faceflux;
        }
    }
    return Fl;
}

template <typename Mesh>
void DarcySolver<Mesh>::postProcessTraceOfPressureRT0()
{
    if(verbose) std::cout << "Postprocessing of TP (RT0 per element)\n";
    if(post_proc_format == "medit"){
        wr_medit_ascii_scalar(post_dir + "/presTP0.bb",globalTP.giveVec(),globalTP.size(),1);
    } else {
        std::cerr
            <<"Warning: Solution constant by element is possible only with medit for the moment\n";
    }
}

template <typename Mesh>
void DarcySolver<Mesh>::postProcessVelocityRT0()
{
    if(verbose) std::cout << "Postprocessing of velocity (RT0 per element)\n";
    if(post_proc_format == "medit"){
        wr_medit_ascii_scalar(post_dir + "/velocRT0.bb",globalFlux.giveVec(),globalFlux.size(),1);
    } else {
        std::cerr
            <<"Warning: Solution constant by element is possible only with medit for the moment\n";
    }
}

template <typename Mesh>
void DarcySolver<Mesh>::postProcessPressureQ0()
{
    if ( verbose > 1 )
        std::cout << "Postprocessing of pressure (constant by element)\n";
    if ( post_proc_format == "medit" ) {
        wr_medit_ascii_scalar(post_dir + "/presQ0.bb",globalP.giveVec(),globalP.size(),1);
    } else {
        std::cerr
            <<"Warning: Solution constant by element is possible only with medit for the moment\n";
    }

#if ANALYTICAL_SOL
    if ( verbose > 1 )
        std::cout <<"Compute L2 pressure error:\n";
    AnalyticalSolPres analyticSol;

    double normL2=0., normL2diff=0., normL2sol=0.;
    double normL2sq=0., normL2diffsq=0., normL2solsq=0.;

    for(UInt i=1; i<=this->_mesh.numVolumes(); ++i){

        pfe.updateFirstDeriv(this->_mesh.volumeList(i));

        normL2sq     += elem_L2_2(globalP,pfe,pdof);
        normL2solsq  += elem_L2_2(analyticSol,pfe);
        normL2diffsq += elem_L2_diff_2(globalP,analyticSol,pfe,pdof);

    }

    normL2     = sqrt(normL2sq);
    normL2sol  = sqrt(normL2solsq);
    normL2diff = sqrt(normL2diffsq);

    std::string errname = post_dir + "/errQ0Pres.txt";
    std::ofstream ofile(errname.c_str());

    ASSERT(ofile,"Error: Output file cannot be opened.");
    ofile << "PRESSION ERROR (Q0)" << std::endl;
    ofile << "|| p       ||_{L^2}                   = " << normL2 << std::endl;
    ofile << "|| p_ex     ||_{L^2}                   = " << normL2sol << std::endl;
    ofile << "|| p - p_ex ||_{L^2}                   = " << normL2diff<< std::endl;
    ofile << "|| P - p_ex ||_{L^2} / || p_ex ||_{L^2} = "
          << normL2diff / normL2sol << "\n" << std::endl;
    ofile << "SQUARE of PRESSION ERROR (Q0)" << std::endl;
    ofile << "|| p       ||^2_{L^2}                   = " << normL2sq << std::endl;
    ofile << "|| p_ex     ||^2_{L^2}                   = " << normL2solsq << std::endl;
    ofile << "|| p - p_ex ||^2_{L^2}                   = " << normL2diffsq << std::endl;
    ofile << "|| P - p_ex ||^2_{L^2} / || p_ex ||^2_{L^2} = "
          << normL2diffsq / normL2solsq << std::endl;
#endif

}

template <typename Mesh>
void DarcySolver<Mesh>::projectPressureQ1( ScalUnknown<Vector> & p_q1 )
{
    // Q1 or P1 elements
    CurrentFE fe_q1(refPFEnodal,geoMap,qr);
    Dof dof_q1(refPFEnodal);
    dof_q1.update(this->_mesh);
    UInt dim_q1 = dof_q1.numTotalDof();
    ScalUnknown<Vector> f_q1(dim_q1);
    p_q1=0.0;
    f_q1=0.0;
    MSRPatt pattA_q1(dof_q1);
    MSRMatr<double> A_q1(pattA_q1);
    ElemMat elmat(fe_q1.nbNode,1,1);
    ElemVec elvec(fe_q1.nbNode,1);
    for(UInt i = 1; i<=this->_mesh.numVolumes(); i++){
        fe_q1.updateJac(this->_mesh.volumeList(i));
        elmat.zero();
        elvec.zero();
        mass(1.,elmat,fe_q1);
        source(globalP(i-1),elvec,fe_q1,0);
        assemb_mat(A_q1,elmat,fe_q1,dof_q1,0,0);
        assemb_vec(f_q1,elvec,fe_q1,dof_q1,0);
    }
    int    options[AZ_OPTIONS_SIZE];
    double params[AZ_PARAMS_SIZE];
    // we first initialize Aztec with its defaults and user's parameters
    aztecOptionsFromDataFile(options,params);
    // next, we overload some of them
    params[AZ_tol] = aztec_tol;
    options[AZ_solver] = AZ_cg;
    options[AZ_precond] = AZ_dom_decomp;
    options[AZ_subdomain_solve] = AZ_icc;
    aztecSolveLinearSyst(A_q1,p_q1.giveVec(),f_q1.giveVec(),p_q1.size(),
                         pattA_q1,options,params);
}

template <typename Mesh>
void DarcySolver<Mesh>::postProcessPressureQ1()
{
    if( verbose > 1 )
        std::cout << "Postprocessing of pressure (L2 projection on the nodes)\n";
    // Q1 or P1 elements
    CurrentFE fe_q1( refPFEnodal , geoMap , qr );
    Dof dof_q1( refPFEnodal );
    dof_q1.update( this->_mesh );
    UInt dim_q1 = dof_q1.numTotalDof();
    ScalUnknown<Vector> nodalPres( dim_q1 );
    projectPressureQ1( nodalPres );

    std::string vtkname,bbname;
    /*
      char str_iter[10],str_time[10];
      static int iter_post=0;
      sprintf(str_time,"t=%f",time);
      sprintf(str_iter,".%03d",iter_post);
    */
    //
    if(post_proc_format == "medit"){
        //bbname = post_dir + "/presQ1" + (string) str_iter + ".bb";
        bbname = post_dir + "/presQ1.bb";
        wr_medit_ascii_scalar(bbname,nodalPres.giveVec(),nodalPres.size());
    } else if (post_proc_format == "vtk"){
        //vtkname = post_dir + "/presQ1" + (string) str_iter + ".vtk";
        vtkname = post_dir + "/presQ1.vtk";
        wr_vtk_ascii_header(vtkname,"Pressure",this->_mesh, dof_q1, fe_q1);
        wr_vtk_ascii_scalar(vtkname,"P",nodalPres.giveVec(), nodalPres.size());
    }
    //iter_post ++;
    //---------------------------------------------

    //! in case of existing analytical solution
#if ANALYTICAL_SOL
    if ( verbose > 1 )
        std::cout <<"Compute the pressure error.\n";

    AnalyticalSolPres analyticSol;

    double normL2=0., normL2diff=0., normL2sol=0.;
    double normH1=0., normH1diff=0., normH1sol=0.;

    for(UInt ivol=1; ivol<=this->_mesh.numVolumes(); ++ivol){

        fe_q1.updateFirstDeriv(this->_mesh.volumeList(ivol));

        normL2     += elem_L2_2(nodalPres,fe_q1,dof_q1);
        normL2sol  += elem_L2_2(analyticSol,fe_q1);
        normL2diff += elem_L2_diff_2(nodalPres,analyticSol,fe_q1,dof_q1);

        normH1     += elem_H1_2(nodalPres,fe_q1,dof_q1);
        normH1sol  += elem_H1_2(analyticSol,fe_q1);
        normH1diff += elem_H1_diff_2(nodalPres,analyticSol,fe_q1,dof_q1);
    }

    normL2     = sqrt(normL2);
    normL2sol  = sqrt(normL2sol);
    normL2diff = sqrt(normL2diff);

    normH1     = sqrt(normH1);
    normH1sol  = sqrt(normH1sol);
    normH1diff = sqrt(normH1diff);

    std::string errname = post_dir + "/errQ1Pres.txt";
    std::ofstream ofile(errname.c_str());

    ASSERT(ofile,"Error: Output file cannot be opened.");
    ofile << "PRESSION ERROR (Q1)" << std::endl;
    ofile << "|| p       ||_{L^2}                   = " << normL2 << std::endl;
    ofile << "|| p_ex     ||_{L^2}                   = " << normL2sol << std::endl;
    ofile << "|| p - p_ex ||_{L^2}                   = " << normL2diff<< std::endl;
    ofile << "|| P - p_ex ||_{L^2} / || p_ex ||_{L^2} = " << normL2diff/normL2sol
          << std::endl;

    ofile << "|| U       ||_{H^1}                   = " << normH1 << std::endl;
    ofile << "|| sol     ||_{H^1}                   = " << normH1sol << std::endl;
    ofile << "|| U - sol ||_{H^1}                   = " << normH1diff<< std::endl;
    ofile << "|| U - sol ||_{H^1} / || sol ||_{H^1} = " << normH1diff/normH1sol
          << std::endl;
#endif
}

template <typename Mesh>
void DarcySolver<Mesh>::projectVelocityQ1( PhysVectUnknown<Vector>& u_q1 )
{
    // Q1 or P1 elements
    CurrentFE fe_q1( refPFEnodal , geoMap , qr );
    Dof dof_q1( refPFEnodal );
    dof_q1.update( this->_mesh );
    UInt dim_q1 = dof_q1.numTotalDof();
    PhysVectUnknown<Vector> f_q1( dim_q1 );
    u_q1=0.0;
    f_q1=0.0;
    MSRPatt pattA_q1(dof_q1,nbCoor);
    MSRMatr<double> A_q1(pattA_q1);
    ElemMat elmat_hdiv(fe_q1.nbNode,nbCoor,0,
                       vfe.nbNode,0,1);
    ElemMat elmat(fe_q1.nbNode,nbCoor,nbCoor);
    ElemVec elvec(fe_q1.nbNode,nbCoor);
    ElemVec elvec_hdiv(vfe.nbNode,1);
    Tab1dView elvec_hdiv_vec = elvec_hdiv.block(0);

    for(UInt i = 1; i<=this->_mesh.numVolumes(); i++){
        fe_q1.updateJac(this->_mesh.volumeList(i));
        vfe.updatePiola(this->_mesh.volumeList(i));
        elmat.zero();
        elmat_hdiv.zero();
        mass(1.,elmat,fe_q1,0,0,nbCoor);
        mass_Mixed_Hdiv(1.,elmat_hdiv,fe_q1,vfe,0,0);
        extract_vec(globalFlux,elvec_hdiv,refVFE,vdof,this->_mesh.volumeList(i).id(),0);
        //
        for(int j=0;j<(int) this->_mesh.volumeList(i).numLocalFaces;j++){
            elvec_hdiv_vec[j] *= signLocalFace( (int)this->_mesh.volumeList(i).id() - 1, j);
        }
        //
        elvec.vec() = elmat_hdiv.mat() * elvec_hdiv.vec();
        for(UInt icoor = 0; icoor < nbCoor; icoor++){
            assemb_mat(A_q1,elmat,fe_q1,dof_q1,icoor,icoor);
            assemb_vec(f_q1,elvec,fe_q1,dof_q1,icoor);
        }
    }
    int    options[AZ_OPTIONS_SIZE];
    double params[AZ_PARAMS_SIZE];
    // we first initialize Aztec with its defaults and user's parameters
    aztecOptionsFromDataFile(options,params);
    // next, we overload some of them
    params[AZ_tol] = aztec_tol;
    options[AZ_solver] = AZ_cg;
    options[AZ_precond] = AZ_dom_decomp;
    options[AZ_subdomain_solve] = AZ_icc;
    aztecSolveLinearSyst(A_q1,u_q1.giveVec(),f_q1.giveVec(),u_q1.size(),
                         pattA_q1,options,params);
}


template <typename Mesh>
void DarcySolver<Mesh>::postProcessVelocityQ1()
{
    if( verbose > 1 )
        std::cout << "Postprocessing of velocity (L2 projection on the nodes)\n";
    // Q1 or P1 elements
    CurrentFE fe_q1( refPFEnodal , geoMap , qr );
    Dof dof_q1( refPFEnodal );
    dof_q1.update( this->_mesh );
    UInt dim_q1 = dof_q1.numTotalDof();
    PhysVectUnknown<Vector> nodalVel( dim_q1 );
    projectVelocityQ1( nodalVel );
    //
    std::string vtkname,bbname;
    /*
      char str_iter[10],str_time[10];
      static int iter_post=0;
      sprintf(str_time,"t=%f",time);
      sprintf(str_iter,".%03d",iter_post);
    */
    //
    if(post_proc_format == "medit"){
        //bbname = post_dir + "/vel" + (string) str_iter + ".bb";
        bbname = post_dir + "/velQ1.bb";
        wr_medit_ascii_vector(bbname,nodalVel.giveVec(),nodalVel.size());
    } else if (post_proc_format == "vtk"){
        //vtkname = post_dir + "/vel" + (string) str_iter + ".vtk";
        vtkname = post_dir + "/velQ1.vtk";
        wr_vtk_ascii_header(vtkname,"Velocity",this->_mesh, dof_q1, fe_q1);
        wr_vtk_ascii_vector(vtkname,"U",nodalVel.giveVec(), nodalVel.size());
    }
    //iter_post ++;

    //---------------------------------------------

#if ANALYTICAL_SOL
    if ( verbose > 1 )
        std::cout <<"Compute the velocity error.\n";
    AnalyticalSolFlux analyticSol;

    /*
      analyticSol.init(diffusion_tensor(0,0)  / diffusion_scalar,
      diffusion_tensor(2,0) / diffusion_scalar,
      diffusion_tensor(2,2) / diffusion_scalar,
      diffusion_tensor(1,1) / diffusion_scalar);
    */

    double normL2=0., normL2diff=0., normL2sol=0.;
    double normH1=0., normH1diff=0., normH1sol=0.;

    for(UInt i=1; i<=this->_mesh.numVolumes(); ++i){

        fe_q1.updateFirstDeriv(this->_mesh.volumeList(i));

        normL2     += elem_L2_2(nodalVel,fe_q1,dof_q1,3);
        normL2sol  = -1.; //! elem_L2_2 is not reckognized
        //! (confusion with another templated function)
        // normL2sol  += elem_L2_2<AnalyticalSolFlux>(analyticSol,fe_q1,0.0,3);
        normL2diff += elem_L2_diff_2(nodalVel,analyticSol,fe_q1,dof_q1,0.,3);
        /*
          normH1     += elem_H1_2(nodalVel,fe_q1,dof_q1,0,3);
          normH1sol  += elem_H1_2(analyticSol,fe_q1,0,3);
          normH1diff += elem_H1_diff_2(nodalVel,analyticSol,fe_q1,dof_q1,0,3);
        */
    }

    normL2     = sqrt(normL2);
    normL2sol  = sqrt(normL2sol);
    normL2diff = sqrt(normL2diff);

    normH1     = sqrt(normH1);
    normH1sol  = sqrt(normH1sol);
    normH1diff = sqrt(normH1diff);

    std::string errname = post_dir + "/errQ1Vel.txt";
    std::ofstream ofile(errname.c_str());

    ASSERT(ofile,"Error: Output file cannot be opened.");
    ofile << "VELOCITY ERROR (Q1)" << std::endl;

    ofile << "|| U         ||_{L^2}                   = " << normL2 << std::endl;
    ofile << "|| exact     ||_{L^2}                   = " << normL2sol << std::endl;
    ofile << "|| U - exact ||_{L^2}                   = " << normL2diff<< std::endl;
    ofile << "|| U - exact ||_{L^2}/|| exact ||_{L^2} = " << normL2diff/normL2sol
          << std::endl;

    //  std::cerr << "|| U       ||_{H^1}                   = " << normH1 << std::endl;
    //  std::cerr << "|| exact     ||_{H^1}                   = " << normH1sol << std::endl;
    //  std::cerr << "|| U - exact ||_{H^1}                   = " << normH1diff<< std::endl;
    // std::cerr << "|| U - exact ||_{H^1} / || exact ||_{H^1} = " << normH1diff/normH1sol
    //   << std::endl;

#endif
}


template <typename Mesh>
void DarcySolver<Mesh>::postProcessEnsight()
{
    if(verbose)
        std::cout << "Postprocessing of pressure and velocity (P1-Output)\n";

    // ********** P1 computation of the velocity **********************
    CurrentFE fe_q1( refPFEnodal , geoMap , qr );
    Dof dof_q1( refPFEnodal );
    dof_q1.update( this->_mesh );
    UInt dim_q1 = dof_q1.numTotalDof();
    ScalUnknown<Vector> nodalPres( dim_q1 );
    projectPressureQ1( nodalPres );
    PhysVectUnknown<Vector> nodalVel( dim_q1 );
    projectVelocityQ1( nodalVel );

    Real time = 0.01; // needed for the index of the result-files
    outensight7Mesh3D( this->_mesh, nodalVel, nodalPres,time );

}

}
#endif
