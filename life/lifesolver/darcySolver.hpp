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

#include <life/lifecore/debug.hpp>

#include <life/lifesolver/darcySolverBase.hpp>

#include <life/lifesolver/darcyHandler.hpp>
#include <life/lifefem/bcManage.hpp>

#include <life/lifealg/clapack.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>
#include <life/lifefilters/ensight7Writer.hpp>
#include <life/lifefem/sobolevNorms.hpp>


namespace LifeV
{


/*!
  \brief A mixed hybrid Darcy solver
  \file darcySolver.hpp
  \author J.-F. Gerbeau and V. Martin
  \date 11/2002

  Templated class to be used both for tetra or hexa meshes.
  can be used as follows:

  \verbatim
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
  \endverbatim
*/
template <typename Mesh>
class DarcySolver
    :
    public DarcyHandler<Mesh>,
    virtual public DarcySolverBase
{
public:

    typedef DarcySolverBase::source_type source_type;
    typedef DarcySolverBase::pressure_solution_type pressure_solution_type;
    typedef DarcySolverBase::velocity_solution_type velocity_solution_type;
    typedef DarcySolverBase::error_signal_type error_signal_type;


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


    /**
       compute the matrix for TP and the contribution from the source
       term to the globalF right hand side.
    */
    void computeHybridMatrixAndSourceRHS();

    void applyBC(); //!< apply the b.c. for the TP problem

    void setup()
        {
            computeHybridMatrixAndSourceRHS();
            applyBC();
        }

     void setBC( BCHandler const& __bch )
        {
            _M_bc = __bch;

            // update the dof with the b.c.
            _M_bc.bdUpdate(this->_mesh, this->feBd, this->tpdof);

            if(this->verbose>2)
                _M_bc.showMe(true);

        }

    void solve();//!< solve the linear system for TP with Aztec
    void computePresFlux();//!< Compute P and U (once TP known)

    //! projection of P (Q0/P0) on Q1/P1.
    void projectPressureQ1( ScalUnknown<Vector>& nodalPres );

    //!< projection of U (RT0) on Q1/P1.
    void projectVelocityQ1( PhysVectUnknown<Vector>& nodalVel );

#if 0
    void postProcessTraceOfPressureRT0();//!< postprocess TP constant per face
    void postProcessVelocityRT0();//!< postprocess Velocity (RT0 per element)
    void postProcessPressureQ0();//!< postprocess P constant per element
    void postProcessPressureQ1(); //!< postproc of Q1/P1 pressure.
    void postProcessVelocityQ1();//!<  postproc of Q1/P1 velocity.
    void postProcessEnsight(); //!< postprocessing in ensight format of P and U
#endif

    Real computeFluxFlag(int flag);

    //! L2 error for pressure wrt analytical solution
    void errorL2( darcy_unknown_type, pressure_solution_type );

    //! L2 error for velocity wrt analytical solution
    void errorL2( velocity_solution_type  );

private:
    void _element_computation(int i); //!< computations of element matrices

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

    //!< boundary conditions handler
    BCHandler _M_bc;

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
    msrPattern(this->tpdof),
    mat(msrPattern),
    globalTP(this->dimTPdof),
    globalF(this->dimTPdof),
    globalP(this->dimPdof),
    globalFlux(this->dimTPdof),
    elvecHyb(this->refTPFE.nbDof,1),
    elvecSource(this->pfe.nbNode,1),
    elvecFlux(this->vfe.nbNode,1),
    elmatMix(this->vfe.nbNode,1,1, this->pfe.nbNode,0,1, this->refTPFE.nbDof,0,1),
    elmatHyb(this->refTPFE.nbDof,1,1),
    BtB(this->pfe.nbNode,this->pfe.nbNode),
    CtC(this->refTPFE.nbDof, this->refTPFE.nbDof),
    BtC(this->pfe.nbNode, this->refTPFE.nbDof),
    signLocalFace( this->_mesh.numVolumes(),this->numFacesPerVolume),
    diffusion_scalar_ele( 1 /*this->_mesh.numVolumes()*/) //to save memory when possible...
{
    signLocalFace = -1.0;
    this->diffusion_scalar_ele = this->diffusion_scalar;

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
    this->vfe.updatePiola(this->_mesh.volumeList(ielem));
    /*
      update the current element for the Pressure (used for the source term).
      The only necessary part is the quadpt.
      The jacobian is not used. (check...)
    */
    this->pfe.updateJacQuadPt(this->_mesh.volumeList(ielem));
    /*
      modify the (0,0) block (A) of the matrix
      the blocks (0,1) (B) and (0,2) (C) are independent of the element
      and have already been computed
    */
    // only one block is set to zero since the other ones are not recomputed
    elmatMix.block(0,0) = 0.;

    // coordinate of the barycenter of the current element
    double xg,yg,zg;
    this->pfe.barycenter(xg,yg,zg);
    mass_Hdiv( this->inversePermeability( xg, yg, zg ), elmatMix,this->vfe,0,0 );

}

template <typename Mesh>
void DarcySolver<Mesh>::computeHybridMatrixAndSourceRHS()
{
    int INFO[1] = {0};
    int NBRHS[1] = {1};//  nb columns of the rhs := 1.
    int NBP[1] = {this->pfe.nbNode};//  pressure dof
    int NBU[1] = {this->vfe.nbNode};//  velocity dof
    int NBL[1] = {this->refTPFE.nbDof};//  trace of pressure dof (lagrange multiplier)
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
    grad_Hdiv(1.,elmatMix,this->vfe,this->pfe,0,1);
    /*
      Update the Boundary Matrix
      (independant of the current element thanks to the Piola transform.)
    */
    TP_VdotN_Hdiv(1., elmatMix, this->refTPFE,this->refVdotNFE, 0, 2 );
    //
    if(this->verbose>3){
        Debug( 6100 ) << "elmatHyb : \n" << "\n";
        elmatHyb.showMe();
        Debug( 6100 ) << "elmatMix : \n" << "\n";
        elmatMix.showMe();
    }
    //
    globalTP=ZeroVector( globalTP.size() );
    globalF=ZeroVector( globalF.size() );
    globalP = ZeroVector( globalP.size() );
    globalFlux = ZeroVector( globalFlux.size() );
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

    ElemMat::matrix_type AA = elmatMix.block(0,0);
    ElemMat::matrix_type BB = elmatMix.block(0,1);
    ElemMat::matrix_type CC = elmatMix.block(0,2);

    Debug( 6100 ) << "number of volumes: " << this->_mesh.numVolumes() << "\n";
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
        ElemVec::vector_view RHSTP = elvecHyb.block(0);

        source( this->sourceTerm(), elvecSource, this->pfe, 0 );

        // initialize the rhs vector (clean this some day...)
        ElemVec::vector_view  rhs = elvecSource.block(0); // corresponds to F2
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
        assemb_mat_symm_lower(mat,elmatHyb,this->refTPFE,this->tpdof,
                              this->_mesh.volumeList(ivol).id(),0,0);
        assemb_vec(globalF,elvecHyb,this->refTPFE,this->tpdof,
                   this->_mesh.volumeList(ivol).id(), 0);
        //-----------------------------------
        // END OF LOOP ON THE VOLUME ELEMENTS
        //-----------------------------------
    }
}

template <typename Mesh>
void DarcySolver<Mesh>::applyBC()
{
    bcManage(mat,globalF,this->_mesh,this->tpdof,_M_bc,this->feBd,1.,0.0);
}

template <typename Mesh>
void DarcySolver<Mesh>::solve()
{
    aztecSolveLinearSyst(mat,globalTP.giveVec(),globalF.giveVec(),
                         globalTP.size(),msrPattern);
    computePresFlux();

    // ********** P1 computation of the velocity **********************
    CurrentFE fe_q1( this->refPFEnodal , this->geoMap , this->qr );
    Dof dof_q1( this->refPFEnodal );
    dof_q1.update( this->_mesh );
    UInt dim_q1 = dof_q1.numTotalDof();
    ScalUnknown<Vector> nodalPres( dim_q1 );
    projectPressureQ1( nodalPres );
    PhysVectUnknown<Vector> nodalVel( dim_q1 );
    projectVelocityQ1( nodalVel );

    //solve_signal_type _M_solve_signal;

    Real time = 0.01; // needed for the index of the result-files
    outensight7Mesh3D( this->_mesh, nodalVel, nodalPres,time );

}

template <typename Mesh>
void DarcySolver<Mesh>::computePresFlux()
{
    //! reset to 0
    globalP = ZeroVector( globalP.size() );
    globalFlux = ZeroVector( globalFlux.size() );

    int INFO[1] = {0};
    int NBRHS[1] = {1};// nb columns of the rhs := 1.
    int INC1[1] = {1};// increment := 1.
    int NBP[1] = {this->pfe.nbNode};//  pressure dof
    int NBU[1] = {this->vfe.nbNode};//  velocity dof
    int NBL[1] = {this->refTPFE.nbDof};//  trace of pressure dof (lagrange multiplier)
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
    ElemMat::matrix_type AA = elmatMix.block(0,0);
    ElemMat::matrix_type BB = elmatMix.block(0,1);
    ElemMat::matrix_type CC = elmatMix.block(0,2);

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
        source( this->sourceTerm(), elvecSource, this->pfe,  0 );

        // initialize the rhs vector (clean this some day...)
        ElemVec::vector_view rhs = elvecSource.block(0); // corresponds to F2
        // Compute rhs = LB^{-1} rhs <- LB^{-1} F2
        dtrtrs_("L", "N", "N", NBP, NBRHS, BtB, NBP, rhs, NBP, INFO);
        ASSERT_PRE(!INFO[0],
                   "Lapack Computation rhs = LB^{-1} rhs is not achieved.");
        // extract the resulting TP for the current fe and put it into elvecHyb.
        extract_vec(globalTP, elvecHyb, this->refTPFE, this->tpdof,ivol, 0);
        ElemVec::vector_view  RHSTP = elvecHyb.block(0);
        // RHSTP = elvecHyb.block(0)  contains the local TP for the current fe.
        // rhs = BtC * RHSTP + rhs <- LB^{-1} Bt A^{-1} C * L + LB^{-1} F2
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
        assemb_vec( globalP, elvecSource, this->refPFE, this->pdof,ivol, 0);
        //__________________________________
        // 2/ Computation of the VELOCITIES
        //__________________________________
        // initialize the element flux vector.
        elvecFlux.zero();     //  Flux (NBU)
        // initialize the flux vector (clean this some day...)
        ElemVec::vector_view flux = elvecFlux.block(0);
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


template <typename Mesh>
void
DarcySolver<Mesh>::errorL2( darcy_unknown_type __type,  pressure_solution_type __analytical_sol )
{
    Debug( 6100 ) <<"Compute L2 pressure error:\n";

    double normL2sq=0.;
    double normL2diffsq=0.;
    double normL2solsq=0.;

    switch( __type )
    {
        case DARCY_PRESSURE_GLOBAL:
        {
            for(UInt i=1; i<=this->_mesh.numVolumes(); ++i)
            {
                this->pfe.updateFirstDeriv(this->_mesh.volumeList(i));


                normL2sq     += elem_L2_2(globalP,this->pfe,this->pdof);
                normL2solsq  += elem_L2_2( __analytical_sol, this->pfe );
                normL2diffsq += elem_L2_diff_2(globalP,__analytical_sol,this->pfe,this->pdof);
            }
        }
        break;
        case DARCY_PRESSURE:
        {

            // Q1 or P1 elements
            CurrentFE fe( this->refPFEnodal , this->geoMap , this->qr );
            Dof dof( this->refPFEnodal );
            dof.update( this->_mesh );
            UInt dim = dof.numTotalDof();
            ScalUnknown<Vector> nodalPres( dim );
            projectPressureQ1( nodalPres );

            Debug( 6100 ) << "Postprocessing of pressure (L2 projection on the nodes)\n";
            for(UInt i=1; i<=this->_mesh.numVolumes(); ++i)
            {

                fe.updateFirstDeriv(this->_mesh.volumeList(i));

                normL2sq     += elem_L2_2( nodalPres, fe, dof );
                normL2solsq  += elem_L2_2( __analytical_sol, fe );
                normL2diffsq += elem_L2_diff_2( nodalPres, __analytical_sol, fe, dof );

            }
        }
        break;
        default:
            std::ostringstream __ex;
            __ex << "invalid darcy unknown type: " << __type;
            throw std::invalid_argument( __ex.str() );
            break;
    }

    // send the signals to all observers that the l2 errors has been computed
    _M_error_signal( "L2", __type, normL2sq, normL2solsq, normL2diffsq );
}

template <typename Mesh>
void
DarcySolver<Mesh>::errorL2( velocity_solution_type __analytical_sol )
{
    Debug( 6100 ) <<"Compute L2 velocity error:\n";
    Debug( 6100 ) << "Postprocessing of velocity (L2 projection on the nodes)\n";

    // Q1 or P1 elements
    CurrentFE fe( this->refPFEnodal , this->geoMap , this->qr );
    Dof dof( this->refPFEnodal );
    dof.update( this->_mesh );
    UInt dim = dof.numTotalDof();
    PhysVectUnknown<Vector> nodalVel( dim );
    projectVelocityQ1( nodalVel );

    double normL2sq=0.;
    double normL2diffsq=0.;
    double normL2solsq=0.;


    for(UInt i=1; i<=this->_mesh.numVolumes(); ++i)
    {
        fe.updateFirstDeriv(this->_mesh.volumeList(i));

        normL2sq     += elem_L2_2( nodalVel, fe, dof, 3 );
        normL2solsq  += elem_L2_2( __analytical_sol, fe, 0.0, 3 );
        normL2diffsq += elem_L2_diff_2( nodalVel, __analytical_sol, fe, dof, 0., 3 );

    }

    // send the signals to all observers that the l2 errors has been computed
    _M_error_signal( "L2", DARCY_VELOCITY, normL2sq, normL2solsq, normL2diffsq );
}
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
void DarcySolver<Mesh>::projectPressureQ1( ScalUnknown<Vector> & p_q1 )
{
    // Q1 or P1 elements
    CurrentFE fe_q1(this->refPFEnodal,this->geoMap,this->qr);
    Dof dof_q1(this->refPFEnodal);
    dof_q1.update(this->_mesh);

    UInt dim_q1 = dof_q1.numTotalDof();

    ScalUnknown<Vector> f_q1(dim_q1);
    p_q1=ZeroVector( dim_q1 );
    f_q1=ZeroVector( dim_q1 );

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
    this->aztecOptionsFromDataFile(options,params);
    // next, we overload some of them
    params[AZ_tol] = this->aztec_tol;
    options[AZ_solver] = AZ_cg;
    options[AZ_precond] = AZ_dom_decomp;
    options[AZ_subdomain_solve] = AZ_icc;
    aztecSolveLinearSyst(A_q1,p_q1.giveVec(),f_q1.giveVec(),p_q1.size(),
                         pattA_q1,options,params);
}

template <typename Mesh>
void DarcySolver<Mesh>::projectVelocityQ1( PhysVectUnknown<Vector>& u_q1 )
{
    // Q1 or P1 elements
    CurrentFE fe_q1( this->refPFEnodal , this->geoMap , this->qr );
    Dof dof_q1( this->refPFEnodal );
    dof_q1.update( this->_mesh );
    UInt dim_q1 = dof_q1.numTotalDof();

    PhysVectUnknown<Vector> f_q1( dim_q1 );
    u_q1=ZeroVector( u_q1.size() ); // warning here number of components is 3 -> size = 3*dim_q1
    f_q1=ZeroVector( f_q1.size() ); // warning here number of components is 3 -> size = 3*dim_q1

    MSRPatt pattA_q1(dof_q1,this->nbCoor);
    MSRMatr<double> A_q1(pattA_q1);
    ElemMat elmat_hdiv(fe_q1.nbNode,this->nbCoor,0,
                       this->vfe.nbNode,0,1);
    ElemMat elmat(fe_q1.nbNode,this->nbCoor,this->nbCoor);
    ElemVec elvec(fe_q1.nbNode,this->nbCoor);
    ElemVec elvec_hdiv(this->vfe.nbNode,1);
    ElemVec::vector_view elvec_hdiv_vec = elvec_hdiv.block(0);

    for(UInt i = 1; i<=this->_mesh.numVolumes(); i++){
        fe_q1.updateJac(this->_mesh.volumeList(i));
        this->vfe.updatePiola(this->_mesh.volumeList(i));
        elmat.zero();
        elmat_hdiv.zero();
        mass(1.,elmat,fe_q1,0,0,this->nbCoor);
        mass_Mixed_Hdiv(1.,elmat_hdiv,fe_q1,this->vfe,0,0);
        extract_vec(globalFlux,elvec_hdiv,this->refVFE,this->vdof,this->_mesh.volumeList(i).id(),0);
        //
        for(int j=0;j<(int) this->_mesh.volumeList(i).numLocalFaces;j++){
            elvec_hdiv_vec[j] *= signLocalFace( (int)this->_mesh.volumeList(i).id() - 1, j);
        }
        //
        elvec.vec() = elmat_hdiv.mat() * elvec_hdiv.vec();
        for(UInt icoor = 0; icoor < this->nbCoor; icoor++){
            assemb_mat(A_q1,elmat,fe_q1,dof_q1,icoor,icoor);
            assemb_vec(f_q1,elvec,fe_q1,dof_q1,icoor);
        }
    }
    int    options[AZ_OPTIONS_SIZE];
    double params[AZ_PARAMS_SIZE];
    // we first initialize Aztec with its defaults and user's parameters
    this->aztecOptionsFromDataFile(options,params);
    // next, we overload some of them
    params[AZ_tol] = this->aztec_tol;
    options[AZ_solver] = AZ_cg;
    options[AZ_precond] = AZ_dom_decomp;
    options[AZ_subdomain_solve] = AZ_icc;
    aztecSolveLinearSyst(A_q1,u_q1.giveVec(),f_q1.giveVec(),u_q1.size(),
                         pattA_q1,options,params);
}

#if 0

//-------------------------------------------------
//! post processing part.
//-------------------------------------------------
//! switch to compute or not the analytical solution
#define ANALYTICAL_SOL 0


template <typename Mesh>
void DarcySolver<Mesh>::postProcessTraceOfPressureRT0()
{
    Debug( 6100 ) << "Postprocessing of TP (RT0 per element)\n";
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
    Debug( 6100 ) << "Postprocessing of velocity (RT0 per element)\n";
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
    Debug( 6100 ) << "Postprocessing of pressure (constant by element)\n";
    if ( post_proc_format == "medit" ) {
        wr_medit_ascii_scalar(post_dir + "/presQ0.bb",globalP.giveVec(),globalP.size(),1);
    } else {
        std::cerr
            <<"Warning: Solution constant by element is possible only with medit for the moment\n";
    }
}


template <typename Mesh>
void DarcySolver<Mesh>::postProcessPressureQ1()
{
    Debug( 6100 ) << "Postprocessing of pressure (L2 projection on the nodes)\n";
    // Q1 or P1 elements
    CurrentFE fe_q1( this->refPFEnodal , this->geoMap , this->qr );
    Dof dof_q1( this->refPFEnodal );
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

}


template <typename Mesh>
void DarcySolver<Mesh>::postProcessVelocityQ1()
{
    Debug( 6100 ) << "Postprocessing of velocity (L2 projection on the nodes)\n";
    // Q1 or P1 elements
    CurrentFE fe_q1( this->refPFEnodal , this->geoMap , this->qr );
    Dof dof_q1( this->refPFEnodal );
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
}


template <typename Mesh>
void DarcySolver<Mesh>::postProcessEnsight()
{
    Debug( 6100 ) << "Postprocessing of pressure and velocity (P1-Output)\n";

    // ********** P1 computation of the velocity **********************
    CurrentFE fe_q1( this->refPFEnodal , this->geoMap , this->qr );
    Dof dof_q1( this->refPFEnodal );
    dof_q1.update( this->_mesh );
    UInt dim_q1 = dof_q1.numTotalDof();
    ScalUnknown<Vector> nodalPres( dim_q1 );
    projectPressureQ1( nodalPres );
    PhysVectUnknown<Vector> nodalVel( dim_q1 );
    projectVelocityQ1( nodalVel );

    Real time = 0.01; // needed for the index of the result-files
    outensight7Mesh3D( this->_mesh, nodalVel, nodalPres,time );

}
#endif
}
#endif
