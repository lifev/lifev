/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s): Miguel A. Fernandez <miguel.fernandez@inria.fr>
Date: 2005-10-09

Copyright (C) 2005 INRIA

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
/*!
  \file NavierStokesAleSolver.h
  \author M.A. Fernandez
  \date 10/2005
  \version 1.0

  \brief This file contains a Navier-Stokes solver allowing moving domains (via an
  ALE formulation). This class implements a time semi-implicit scheme with stabilized
  finite elements in space (SDFEM and IP). Tested with P1/P1 (SDFEM and IP) and Q1/Q1 (only SDFEM)

  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  \note For SDFEM we assume the source to be zero (to be changed...)
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

*/
#ifndef _NAVIERSTOKESALESOLVER_H_
#define _NAVIERSTOKESALESOLVER_H_

#include <life/lifesolver/NavierStokesHandler.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifefem/assembGeneric.hpp>
#include <life/lifefem/bcManage.hpp>

//#include <life/lifefem/meshMotion.hpp>

//#include <life/lifealg/SolverAztec.hpp>

#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifefem/sobolevNorms.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifesolver/sdStabilization.hpp>
#include <life/lifesolver/ipStabilization.hpp>

namespace LifeV
{
    /*!
      \class NavierStokesAleSolver

      This class implements a time semi-implicit scheme with stabilized
      finite elements in space (SDFEM and IP). Tested with P1/P1 (SDFEM and IP) and Q1/Q1 (only SDFEM)
    */

    template<typename Mesh>
    class NavierStokesAleSolver:
        public NavierStokesSolver<Mesh>
    {

    public:
        typedef  typename  NavierStokesHandler<Mesh>::Function Function;
        typedef  typename  NavierStokesHandler<Mesh>::source_type source_type;

        //! Constructor
        /*!
          \param dataFile GetPot data file
          \param refFE reference FE
          \param quadRule volumic quadrature rule
          \param boundaryQuadRule surface quadrature rule
          \param bcHandler boundary conditions for the velocity
        */
        NavierStokesAleSolver( const GetPot& dataFile,
                               const RefFE& refFE,
                               const QuadRule& quadRule,
                               const QuadRule& boundaryQuadRule,
                               BCHandler& bcHandler,
                               BCHandler& bcExtOper);

        NavierStokesAleSolver( const GetPot& dataFile,
                               const RefFE& refFE,
                               const QuadRule& quadRule,
                               const QuadRule& boundaryQuadRule);
        //! Update the right hand side for time advancing
        /*!
          \param source volumic source
          \param time present time
        */
        void timeAdvance( source_type const& source, Real const& time );

        //! Update convective term, bc treatment and solve the linearized ns system
        void iterate( const Real& time );
        void iterateLin( const Real& time, BCHandler& bcHdu, const bool shapeTerms=1);
        void iterateReducedLin( const Real& time, BCHandler& bcHdp, BCHandler& bcHdz );

        //! removes mean of component comp of vector x
        void removeMean( Vector& x, UInt comp=1 );

        //!
//        void updateMesh( const Real& time );
//        void updateDispVelo();

//       void setFullEssential(bool _full)
//           {
//               _factor_data.setFullEssential(_full);
//               _factor_data_jacobian.setFullEssential(_full);
//           }

//         Vector& w() {return M_w;}
//         Vector& dw() {return M_dw;}

//        const Vector& d() const {return M_d;}

        Vector& residual() {return M_res;}
        Vector  getDeltaLambda() {return this->dt()*M_dsol;}

//         Vector& dp() {return M_dp;}

    private:

        //! constructor setup
        void     setUp(const GetPot&);

        //! Block pattern of full matrix
        MSRPatt M_fullPattern, M_uPattern;

        //! Matrix M_u: Vmass
        //MixedMatr<nDimensions, nDimensions, MSRPatt, double> M_matrMass;
        MSRMatr<double> M_matrMass;

        //! Matrix C: CStokes + Convective_term + nonlinear stabilizations
        MSRMatr<double> M_matrFull, M_matrFullAux;

        //! Elementary matrices and vectors
        ElemMat M_elmatC; //velocity stiffness
        ElemMat M_elmatMass; //velocity mass
        ElemMat M_elmatD;
        ElemMat M_elmatDtr;
        ElemVec M_elvec; // Elementary right hand side
        ElemVec M_uLoc, M_wLoc; // Elementary mesh velocity

        Vector M_beta;

        //! Right  hand side for the velocity
        Vector M_rhsU;

        //! Right  hand side global
        Vector M_rhsFull, M_rhsFullAux;



        //! Global solution _u and _p
        Vector M_sol, M_res;

        //! Backward velocity and mesh velocity
        Vector M_un;
//        Vector  M_dOld;

        SolverAztec M_linearSolver;

        Real M_time;

        Vector M_constantPressure;

//        HarmonicExtension M_extOper;
//        const Vector& M_d;

        SDStabilization<Mesh, Dof> M_sdStab;

        IPStabilization<Mesh, Dof> M_ipStab;

        // terms for shape derivative

        ElemVec M_elvec_dudp; // Elementary right hand side for the linearized velocity
        ElemVec M_w_loc; // Elementary mesh velocity
        ElemVec M_uk_loc; // Elementary velocity
        ElemVec M_pk_loc; // Elementary pressure
        ElemVec M_d_loc; // Elementary displacement for right hand side
        ElemVec M_dw_loc; // Elementary mesh velocity for right hand side

        Vector  M_rhsFullAuxdz;
        Vector  M_dsol;
//        Vector  M_dw;

        void M_computeVelocityMassMatrix();

        bool M_firstTimeStep;


        // Reduced fluid solver for Newton FSI
        MSRPatt M_pPattern;
        MSRMatr<double> M_matrL, M_matrLAux;
        bool M_newDomain;
        SolverAztec M_linearSolverReduced;
        Vector M_dp;
        UInt M_nUsePC;
        Real M_tTotalSolve;
        Real M_tLastSolve;
        Real M_tThisSolve;

        /*! solve linear system
          @return number of iterations, zero if not applicable
          @param condEst condition estimate is returned here, -1 if not applicable
        */
        UInt solveLinearSystem( Real& condEst );

        //! solve linear system once (iterative solvers)
        void solveLinearSystemOnce( bool reusePC );

    }; // class NavierStokesAleSolver



    //
    //                                         IMPLEMENTATION
    //
    template<typename Mesh> NavierStokesAleSolver<Mesh>::
    NavierStokesAleSolver( const GetPot& dataFile,
                           const RefFE& refFE,
                           const QuadRule& quadRule,
                           const QuadRule& boundaryQuadRule,
                           BCHandler& bcHandler,
                           BCHandler& bcExtOper ):
        NavierStokesAleHandler<Mesh>( dataFile, refFE, refFE, quadRule,
                                   boundaryQuadRule, quadRule, boundaryQuadRule,
                                   bcHandler ),
        M_fullPattern( this->_dof_u, this->_mesh, nDimensions+1, this->patternType() ),
        M_uPattern( this->_dof_u, this->_mesh, nDimensions, this->patternType() ),
        M_matrMass( M_uPattern ),
        M_matrFull( M_fullPattern ),
        M_matrFullAux( M_fullPattern ),
        M_elmatC( this->_fe_u.nbNode, nDimensions, nDimensions ),
        M_elmatMass( this->_fe_u.nbNode, nDimensions, nDimensions ),
        M_elmatD( this->_fe_u.nbNode, nDimensions+1, nDimensions ),
        M_elmatDtr( this->_fe_u.nbNode, nDimensions, nDimensions+1 ),
        M_elvec( this->_fe_u.nbNode, nDimensions+1 ),
        M_uLoc( this->_fe_u.nbNode, nDimensions ),
        M_wLoc( this->_fe_u.nbNode, nDimensions ),
        M_beta( nDimensions*this->_dim_u ),
        M_rhsU( nDimensions*this->_dim_u ),
        M_rhsFull( ( nDimensions+1 )*this->_dim_u ),
        M_rhsFullAux( ( nDimensions+1 )*this->_dim_u ),
        M_sol( ( nDimensions+1 )*this->_dim_u ),
        M_res( ( nDimensions+1 )*this->_dim_u ),
        M_un( ( nDimensions )*this->_dim_u ),
//        M_w(  ( nDimensions )*this->_dim_u ),
//        M_dOld(  ( nDimensions )*this->_dim_u ),
        M_constantPressure( ( nDimensions+1 )*this->_dim_u ),
//         M_extOper(  dataFile,
//                     this->mesh(), 1.0,
//                     quadRule,
//                     boundaryQuadRule,
//                     bcExtOper ),
//        M_d(  M_extOper.getDisplacement()  ),
        M_sdStab( dataFile, this->mesh(), this->_dof_u, this->_refFE_u, this->_Qr_u, this->viscosity() ),
        M_ipStab( dataFile, this->mesh(), this->_dof_u, this->_refFE_u, this->_feBd_u, this->_Qr_u, this->viscosity() ),
        M_elvec_dudp( this->_fe_u.nbNode, nDimensions + 1),
        M_w_loc( this->_fe_u.nbNode, nDimensions ),
        M_uk_loc( this->_fe_u.nbNode, nDimensions ),
        M_pk_loc( this->_fe_u.nbNode, 1 ),
        M_d_loc( this->_fe_u.nbNode, nDimensions ),
        M_dw_loc( this->_fe_u.nbNode, nDimensions ),
        M_rhsFullAuxdz( ( nDimensions+1 )*this->_dim_u ),
        M_dsol( (nDimensions+1)*this->_dim_u ),
//        M_dw( nDimensions*this->_dim_u ),
        M_firstTimeStep(1),
        M_pPattern( this->_dof_u, this->mesh(), 1, BasePattern::STANDARD_PATTERN ),
        M_matrL( M_pPattern ),
        M_matrLAux( M_pPattern ),
        M_newDomain( 0 ),
        M_dp(this->_dim_u ),
        M_nUsePC( 1 ),
        M_tTotalSolve( 0 ),
        M_tLastSolve( 1 ),
        M_tThisSolve( 1 )
    {

        setUp(dataFile);

    } // NavierStokesAleSolver() [Constructor]


    template<typename Mesh> NavierStokesAleSolver<Mesh>::
    NavierStokesAleSolver( const GetPot&   dataFile,
                           const RefFE&    refFE,
                           const QuadRule& quadRule,
                           const QuadRule& boundaryQuadRule):
        NavierStokesAleHandler<Mesh>( dataFile,
                                      refFE,
                                      refFE,
                                      quadRule,
                                      boundaryQuadRule,
                                      quadRule,
                                      boundaryQuadRule),
        M_fullPattern( this->_dof_u, this->mesh(), nDimensions+1, this->patternType() ),
        M_uPattern( this->_dof_u, this->mesh(), nDimensions, this->patternType() ),
        M_matrMass( M_uPattern ),
        M_matrFull( M_fullPattern ),
        M_matrFullAux( M_fullPattern ),
        M_elmatC( this->_fe_u.nbNode, nDimensions, nDimensions ),
        M_elmatMass( this->_fe_u.nbNode, nDimensions, nDimensions ),
        M_elmatD( this->_fe_u.nbNode, nDimensions+1, nDimensions ),
        M_elmatDtr( this->_fe_u.nbNode, nDimensions, nDimensions+1 ),
        M_elvec( this->_fe_u.nbNode, nDimensions+1 ),
        M_uLoc( this->_fe_u.nbNode, nDimensions ),
        M_wLoc( this->_fe_u.nbNode, nDimensions ),
        M_beta( nDimensions*this->_dim_u ),
        M_rhsU( nDimensions*this->_dim_u ),
        M_rhsFull( ( nDimensions+1 )*this->_dim_u ),
        M_rhsFullAux( ( nDimensions+1 )*this->_dim_u ),
        M_sol( ( nDimensions+1 )*this->_dim_u ),
        M_res( ( nDimensions+1 )*this->_dim_u ),
        M_un( ( nDimensions )*this->_dim_u ),
//        M_w(  ( nDimensions )*this->_dim_u ),
//        M_dOld(  ( nDimensions )*this->_dim_u ),
        M_constantPressure( ( nDimensions+1 )*this->_dim_u ),
//         M_extOper(  dataFile,
//                     this->mesh(), 1.0,
//                     quadRule,
//                     boundaryQuadRule,
//                     bcExtOper ),
//        M_d(  M_extOper.getDisplacement()  ),
        M_sdStab( dataFile, this->mesh(), this->_dof_u, this->_refFE_u, this->_Qr_u, this->viscosity() ),
        M_ipStab( dataFile, this->mesh(), this->_dof_u, this->_refFE_u, this->_feBd_u, this->_Qr_u, this->viscosity() ),
        M_elvec_dudp( this->_fe_u.nbNode, nDimensions + 1),
        M_w_loc( this->_fe_u.nbNode, nDimensions ),
        M_uk_loc( this->_fe_u.nbNode, nDimensions ),
        M_pk_loc( this->_fe_u.nbNode, 1 ),
        M_d_loc( this->_fe_u.nbNode, nDimensions ),
        M_dw_loc( this->_fe_u.nbNode, nDimensions ),
        M_rhsFullAuxdz( ( nDimensions+1 )*this->_dim_u ),
        M_dsol( (nDimensions+1)*this->_dim_u ),
//        M_dw( nDimensions*this->_dim_u ),
        M_firstTimeStep(1),
        M_pPattern( this->_dof_u, this->mesh(), 1, BasePattern::STANDARD_PATTERN ),
        M_matrL( M_pPattern ),
        M_matrLAux( M_pPattern ),
        M_newDomain( 0 ),
        M_dp(this->_dim_u ),
        M_nUsePC( 1 ),
        M_tTotalSolve( 0 ),
        M_tLastSolve( 1 ),
        M_tThisSolve( 1 )
    {
        setUp(dataFile);
    }


    template<typename Mesh> void NavierStokesAleSolver<Mesh>::
    setUp(const GetPot& dataFile)
    {

        // check mesh for elements with all nodes on the boundary
        UInt nLocalFaces = this->mesh().numLocalFaces();
        UInt nFixedTets = 0;
        for ( UInt iVol = 1; iVol <= this->mesh().numVolumes(); iVol++ )
        {
            UInt nBoundaryFaces = 0;
            for ( UInt iLocalFace=1; iLocalFace<=nLocalFaces; ++iLocalFace )
            {
                if ( this->mesh().isBoundaryFace( this->mesh().localFaceId( iVol, iLocalFace )))
                {
                    ++nBoundaryFaces;
                }
            }
            if ( nBoundaryFaces > 1 )
            {
                ++nFixedTets;
            }
        }
        if ( nFixedTets > 0 )
        {
            std::cerr << "WARNING: " << nFixedTets << " of " << this->mesh().numVolumes()
                      << " tetrahedrons have all nodes on the boundary."
                      << std::endl;
        }


        M_linearSolver.setOptionsFromGetPot( dataFile, "fluid/aztec" );
        M_linearSolverReduced.setOptionsFromGetPot( dataFile, "fluid/reduced_aztec" );

        std::cout << std::endl;
        std::cout << "O-  Pressure unknowns: " << this->_dim_p     << std::endl;
        std::cout << "F-  Velocity unknowns: " << this->_dim_u     << std::endl
                  << std::endl;


        M_sol = ZeroVector( M_sol.size() );
        M_res = ZeroVector( M_res.size() );
        this->_p = ZeroVector( this->_p.size()  );
    }


    template<typename Mesh> void NavierStokesAleSolver<Mesh>::
    timeAdvance( source_type const& source, Real const& time )
    {

        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // source parameter has to be removed (moving domains) to be done
        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


        M_time = time;

        std::cout << std::endl;
        std::cout << "F== Now we are at time "<< M_time << " s." << std::endl;

        // Number of velocity components
        UInt nbCompU = this->_u.nbcomp();

        Chrono chrono;
        chrono.start();

        // compute mass matrix at the first time step
        if ( M_firstTimeStep )
        {
            std::cout << "  F-  Computing mass matrix and updating mass and source term on right hand side... "
                      << std::flush;
            M_computeVelocityMassMatrix();
            M_firstTimeStep = 0;
        }
        else
            std::cout << "  F-  Updating mass term on right hand side... " << std::flush;

        M_rhsU = M_matrMass * this->_u;

        chrono.stop();
        std::cout << "done in " << chrono.diff() << " s." << std::endl;

        // Full right hand-side
        for ( UInt iDof = 0; iDof<nDimensions*this->_dim_u; ++iDof )
            M_rhsFullAux[ iDof ] =  M_rhsU[ iDof ];
        for ( UInt iDof = 0; iDof<this->_dim_u; ++iDof )
            M_rhsFullAux[ iDof+nDimensions*this->_dim_u ]= 0.0;

        M_un = this->_u;
        this->_dispOld = this->harmonicExtension().getDisplacement();

    } // timeAdvance()


    template<typename Mesh>
    void NavierStokesAleSolver<Mesh>::iterate( const Real& time )
    {

        //
        //
        //

        M_beta = this->density() * ( M_un - this->_w );

        M_rhsFull = M_rhsFullAux;

        std::cout << "F-  Updating matrices... " << std::flush;

        Chrono chrono;
        chrono.start();

        // Matrices initialization
        M_matrMass.zeros();
        M_matrFullAux.zeros();

        // Number of velocity components
        UInt nbCompU = this->_u.nbcomp();

        UInt iloc, ig;

        Real rhodti=  this->density()/this->dt();

        // Elementary computation and matrix assembling
        // Loop on elements
        for ( UInt iVol = 1; iVol <= this->mesh().numVolumes(); iVol++ )
        {
            this->_fe_u.updateFirstDeriv( this->mesh().volumeList( iVol ) );

            M_elmatC.zero();
            M_elmatMass.zero();
            M_elmatD.zero();
            M_elmatDtr.zero();

            // viscous term
            stiff_strain( 2.0*this->viscosity(), M_elmatC, this->_fe_u );

            // mass term
            mass( rhodti, M_elmatMass, this->_fe_u, 0, 0, nDimensions );
            M_elmatC.mat() += M_elmatMass.mat();

            // convective term
            UInt eleID = this->_fe_u.currentId();
            // Non linear term, Semi-implicit approach
            // M_elvec contains the velocity values in the nodes
            for ( UInt iNode = 0 ; iNode<( UInt )this->_fe_u.nbNode ; iNode++ )
            {
                iloc = this->_fe_u.patternFirst( iNode );
                for ( UInt iComp = 0; iComp<nbCompU; ++iComp )
                {
                    ig = this->_dof_u.localToGlobal( eleID, iloc+1 )-1+iComp*this->_dim_u;
                    M_elvec.vec()[ iloc+iComp*this->_fe_u.nbNode ] = M_beta[ig];
                    M_uLoc.vec() [ iloc + iComp * this->_fe_u.nbNode ] = M_un[ig];
                    M_wLoc.vec() [ iloc + iComp * this->_fe_u.nbNode ] = this->_w[ig];
                }
            }

            // ALE term: - rho div w u v
            mass_divw( -this->density(), M_wLoc,  M_elmatC , this->_fe_u, 0, 0, nbCompU );

            // ALE stab implicit: 0.5 rho div w u v
            mass_divw( 0.5*this->density(), M_uLoc,  M_elmatC , this->_fe_u, 0, 0, nbCompU );

            // compute local convective term
            for ( UInt iComp = 0;iComp<nbCompU;iComp++ )
                for  ( UInt jComp = 0;jComp<nbCompU;jComp++ )
                    grad( jComp, M_elvec, M_elmatC, this->_fe_u, this->_fe_u, iComp, iComp );

            // Assembling
            for ( UInt iComp = 0;iComp<nbCompU;iComp++ )
            {
                // mass + convection + stiffness
                for ( UInt jComp = 0;jComp<nbCompU;jComp++ )
                    assemb_mat( M_matrFullAux, M_elmatC, this->_fe_u, this->_dof_u, iComp,jComp );

                // mass term
                assemb_mat( M_matrMass, M_elmatMass, this->_fe_u, this->_dof_u, iComp, iComp);

                // - ( p , \div v)
                grad( iComp, 1.0, M_elmatDtr, this->_fe_u, this->_fe_u, iComp, nDimensions ); // minus inside

                //   ( q , \div u)
                div(  iComp, -1.0, M_elmatD  , this->_fe_u, this->_fe_u, nbCompU, iComp ); // minus inside

                // assembling
                assemb_mat( M_matrFullAux, M_elmatDtr, this->_fe_u, this->_dof_u, iComp, nbCompU );
                assemb_mat( M_matrFullAux, M_elmatD, this->_fe_u, this->_dof_u, nbCompU, iComp );

            }
        }


        chrono.stop();

        std::cout << "done in " << chrono.diff() << " s."
                  << std::endl;

        chrono.start();

        switch ( this->stabilization() )
        {
            case IP_STABILIZATION:
                std::cout << "  o-  Updating IP stabilization terms (matrix)...          "
                          << std::flush;
                M_ipStab.apply( M_matrFullAux, M_beta );
                break;
            case SD_STABILIZATION:
                std::cout << "  o-  Updating SD stabilization terms (matrix)...          "
                          << std::flush;
                M_sdStab.apply(this->dt(), M_matrFullAux, M_beta );

                // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                // stabilization source term here, to be done...
                // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                //
                //	std::cout << "  o-  Updating SD stabilization terms (rhs)...          "
                //	  << std::flush;
                //M_sdStab.apply(this->_dt,M_rhsFullAux, M_beta, source, M_time);
                break;
            default:
                ERROR_MSG("This stabilization is not yet implemented");
        }

        chrono.stop();
        std::cout << "done in " << chrono.diff() << " s." << std::endl;

        // for BC treatment ( done at each time-step )
        std::cout << "  o-  Applying boundary conditions...          "
                  << std::flush;
        chrono.start();


        M_matrFull = M_matrFullAux;

        // BC manage for the velocity
        this->bcHandler().bdUpdate( this->mesh(), this->_feBd_u, this->_dof_u );
        bcManage( M_matrFull, M_rhsFull, this->mesh(), this->_dof_u,  this->bcHandler(),
                  this->_feBd_u, 1.0, M_time );


        std::cout << "  norm_inf (_f_du) after BC= " << norm_2( M_rhsFull ) << std::endl;
        std::cout << "  norm_inf ( difference ) after BC= " << norm_2( M_rhsFull - M_rhsFullAux ) << std::endl;

        if ( this->bcHandler().hasOnlyEssential() )
            M_matrFull.diagonalize( nDimensions*this->_dim_u, 1.0, M_rhsFull, 0);

        chrono.stop();
        std::cout << "done in " << chrono.diff() << " s." << std::endl;

        // ---------------
        // Solving linear system
        // ---------------

        M_sol = ZeroVector( M_sol.size() );


        std::cout << "  o-  Solving system...                        "
                  << std::flush;
        chrono.start();
        Real condEst;
        UInt iter = solveLinearSystem( condEst );
        chrono.stop();
        std::cout << "done in " << chrono.diff() << " s." << std::endl;

        if ( iter > 0 )
            std::cout << "        number of iterations                        = "
                      << iter << std::endl;
        if ( condEst >= 0 )
            std::cout << "        estimated condition number (preconditioned) = "
                      << condEst << std::endl;

        for ( UInt iDof = 0; iDof<nDimensions*this->_dim_u; ++iDof )
            this->_u[ iDof ] = M_sol[ iDof ];

        for ( UInt iDof = 0; iDof<this->_dim_u; ++iDof )
            this->_p[ iDof ] = M_sol[ iDof+nDimensions*this->_dim_u ];

        M_res = M_rhsFullAux - M_matrFullAux * M_sol;
        std::cout << "norm( M_res )  = " << norm_inf( M_res ) << std::endl;

    } // iterate()

    template <typename Mesh>
    void NavierStokesAleSolver<Mesh>::removeMean( Vector& x, UInt comp )
    {
        Real sum1 = 0.;
        Real sum0 = 0.;
        for ( UInt iVol = 1; iVol <= this->mesh().numVolumes(); iVol++ )
        {
            this->_fe_p.updateFirstDeriv( this->mesh().volumeList( iVol ) );
            sum1 += elem_integral( x, this->_fe_p, this->_dof_p, comp );
            sum0 += this->_fe_p.measure();
        }
        Real mean = sum1/sum0;
        for ( UInt iNode = (comp-1)*this->_dim_u; iNode<comp*this->_dim_u; ++iNode )
        {
            x[ iNode ] -= mean;
        }

    } // removeMean()

//     template <typename Mesh>
//     void NavierStokesAleSolver<Mesh>::updateMesh( const Real& time )
//     {

//         // Updating mesh displacement and velocity
//         this->M_extOper.updateExtension( this->mesh(), time );

//         // Updating mesh points
//         this->mesh().moveMesh( M_d );

//         Real dti = 1.0 / this->_dt;
//         M_w = ( M_d - M_dOld ) * dti;

//         std::cout << "norm( M_d )  = " << norm_inf( M_d ) << std::endl;
//         std::cout << "norm( M_w    )  = " << norm_inf( M_w ) << std::endl;

//         M_newDomain = 1;
//     }

    //
    // Linearized iterate for Newton FSI
    //
    template <typename Mesh>
    void NavierStokesAleSolver<Mesh>::
    iterateLin( const Real& time, BCHandler& bcHdu, const bool shapeTerms)
    {

        Chrono chrono;
        chrono.start();




        //initialize right hand side
        M_rhsFullAuxdz = ZeroVector(  M_rhsFullAuxdz.size() );

        // Number of velocity components
        UInt nbCompU = this->_u.nbcomp();
        UInt iloc, ig;

        std::cout << norm_inf(M_un);
        std::cout << norm_inf(_w);
        std::cout << norm_inf(_dw);
        std::cout << norm_inf(_u);
        std::cout << norm_inf(this->harmonicExtension().getDisplacement()); // d local


        if (shapeTerms) {
            /*

            @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            No taking into account source terms
            i.e f = 0. To be modified
            @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            */

            std::cout << "    F-LINEAR-  Updating right hand side (shape derivative terms)... ";
            UInt iloc,ig,icloc;

            chrono.start();

            // Elementary computation and right hand-side assembling
            // Loop on elements
            for ( UInt i = 1; i <= this->mesh().numVolumes(); ++i )
            {

                this->_fe_u.updateFirstDerivQuadPt( this->mesh().volumeList( i ) );

                // initialization of elementary vectors
                M_elvec_dudp.zero();
                for ( UInt k = 0 ; k < ( UInt ) this->_fe_u.nbNode ; ++k )
                {
                    iloc = this->_fe_u.patternFirst( k );
                    for ( UInt ic = 0; ic < nbCompU; ++ic )
                    {
                        ig = this->_dof_u.localToGlobal( i, iloc + 1 ) - 1 + ic * this->_dim_u;
                        icloc = iloc + ic * this->_fe_u.nbNode ;
                         M_elvec.vec() [ icloc  ] = M_un( ig ) - this->_w( ig );  // u^n - w^k local
                         M_w_loc.vec() [ icloc  ] = this->_w( ig );               // w^k local
                         M_uk_loc.vec()[ icloc  ] = this->_u( ig );          // u^k local
                         M_d_loc.vec() [ icloc  ] = this->harmonicExtension().getDisplacement()( ig ); // d local
                         M_dw_loc.vec()[ icloc  ] = this->_dw( ig );              // dw local
                         M_uLoc.vec()  [ icloc  ] = M_un( ig );              // u^n local
                    }
                    ig = this->_dof_u.localToGlobal( i, iloc + 1 ) - 1;
                    M_pk_loc[ iloc ] = this->_p( ig );  // p^k local
                }

                //
                // Elementary vectors
                //
                //  - \rho ( -\grad w^k :[I\div d - (\grad d)^T] u^k + ( u^n-w^k )^T[I\div d - (\grad d)^T] (\grad u^k)^T , v  )
                source_mass1( -this->density(), M_uk_loc,  M_elvec, M_d_loc, M_elvec_dudp, this->_fe_u );

                //  + \rho * ( \grad u^k dw, v  )
                source_mass2( this->density(), M_uk_loc, M_dw_loc, M_elvec_dudp, this->_fe_u );

                //  - \rho/2 ( \nabla u^n:[2 * I\div d - (\grad d)^T]  u^k , v  )
                source_mass3( -0.5*this->density(), M_uLoc, M_uk_loc, M_d_loc, M_elvec_dudp, this->_fe_u );


                //  - ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  )
                source_stress( -1.0, this->viscosity(), M_uk_loc, M_pk_loc, M_d_loc, M_elvec_dudp, this->_fe_u, this->_fe_u );

                // + \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )
                source_stress2( this->viscosity(), M_uk_loc, M_d_loc, M_elvec_dudp, this->_fe_u );

                //  + ( (\grad u^k):[I\div d - (\grad d)^T] , q  )
                source_press( -1.0, M_uk_loc, M_d_loc, M_elvec_dudp, this->_fe_u, this->_fe_u, nDimensions );
                //
                // Assembling
                //

                // assembling presssure right hand side
                assemb_vec( M_rhsFullAuxdz, M_elvec_dudp, this->_fe_u, this->_dof_u, nDimensions );

                // loop on velocity components
                for ( UInt ic = 0; ic < nbCompU; ic++ )
                {
                    // assembling velocity right hand side
                    assemb_vec( M_rhsFullAuxdz, M_elvec_dudp, this->_fe_u, this->_dof_u, ic );
                }
            }

            chrono.stop();
            std::cout << "done in " << chrono.diff() << "s." << std::endl;
        }

        std::cout << "norm( M_rhsFullAuxdz )  = " << norm_inf( M_rhsFullAuxdz ) << std::endl;

        std::cout << "    F-LINEAR-  Applying boundary conditions... ";
        chrono.start();

        M_matrFull  =  M_matrFullAux;
        M_rhsFull = M_rhsFullAuxdz;

        bcHdu.bdUpdate( this->mesh(), this->_feBd_u, this->_dof_u );
        bcManage( M_matrFull, M_rhsFull, this->mesh(), this->_dof_u,  bcHdu , this->_feBd_u, 1.0, time );
        chrono.stop();
        std::cout << "done in " << chrono.diff() << " s." << std::endl;

        std::cout << "  norm_inf (_f_du) after BC= " << norm_inf( M_rhsFull ) << std::endl;
        std::cout << "  norm_inf ( difference ) after BC= " << norm_inf( M_rhsFull - M_rhsFullAuxdz ) << std::endl;


        // ---------------
        // C * V = F
        // ---------------

        // Zero as initial guess
        M_dsol = ZeroVector(  M_dsol.size() );



        std::cout << "  o-  Solving system... "
                  << std::flush;
        chrono.start();
        M_linearSolver.setRecursionLevel( 1 );

        M_linearSolver.solve( M_dsol, M_rhsFull, SolverAztec::SAME_PRECONDITIONER );


        chrono.stop();
        std::cout << "done in " << chrono.diff() << " s." << std::endl;

        if ( !M_linearSolver.converged() )
        {
            std::cerr << "        WARNING: Solver failed to converge."
                      << std::endl;
        }
        std::cout << "        number of iterations                        = "
                  << M_linearSolver.iterations() << std::endl;

        M_res = M_rhsFullAuxdz - M_matrFullAux * M_dsol;
        std::cout << "norm( M_res_linear )  = " << norm_inf( M_res ) << std::endl;
    }



    //  Updating variations for grid velocity and displacement
//     template <typename Mesh>
//     void NavierStokesAleSolver<Mesh>::
//     updateDispVelo()
//     {
//         // Updating mesh displacement and velocity
//         this->M_extOper.updateExtension( this->mesh(), M_time, 1 );

//         Real dti = 1.0 / this->_dt;

//         std::cout << " max norm d^f = " << norm_inf( M_d ) << std::endl;

//         M_dw = M_d * dti;

//         std::cout << " max norm dw = " << norm_inf( M_dw ) << std::endl;
//     }


    // Compute the velocity mas Matrix at the initialization
    template <typename Mesh> void NavierStokesAleSolver<Mesh>::
    M_computeVelocityMassMatrix()
    {

        // Matrices initialization
        M_matrMass.zeros();

        // Number of velocity components
        UInt nbCompU = this->_u.nbcomp();

        Real rhodti = this->density()/this->dt();

        // Elementary computation and matrix assembling
        // Loop on elements
        for ( UInt iVol = 1; iVol <= this->mesh().numVolumes(); iVol++ )
        {
            this->_fe_u.updateJacQuadPt( this->mesh().volumeList( iVol ) );
            M_elmatMass.zero();

            // mass
            mass(rhodti, M_elmatMass, this->_fe_u, 0, 0, nDimensions );

            for ( UInt iComp = 0;iComp<nbCompU;iComp++ )
                assemb_mat( M_matrMass, M_elmatMass, this->_fe_u, this->_dof_u, iComp, iComp);

        }
    }

    //
    // Reduced Linearized iterate for Newton FSI
    //
    template <typename Mesh>
    void NavierStokesAleSolver<Mesh>::
    iterateReducedLin( const Real& time, BCHandler& bcHdp, BCHandler& bcHdz )
    {

        Chrono chrono;
        chrono.start();


        chrono.start();
        if (M_newDomain) { // computing laplacien fe matrix (for "non-linear" fluid evaluation)

            std::cout << "    F-REDUCED_LINEAR-  Updating matrix... ";

            // Initializing matrix
            M_matrLAux.zeros();

            // Loop on elements
            for ( UInt i = 1; i <= this->mesh().numVolumes(); i++ )
            {

                this->_fe_u.updateFirstDerivQuadPt( this->mesh().volumeList( i ) );

                M_elmatC.zero();

                stiff(1.0, M_elmatC, this->_fe_u );
                assemb_mat( M_matrLAux, M_elmatC, this->_fe_u, this->_dof_u, 0, 0);

            }
        }

        chrono.stop();
        std::cout << "done in " << chrono.diff() << "s." << std::endl;

        std::cout << "    F-REDUCED_LINEAR-  Applying boundary conditions... ";
        chrono.start();

        M_matrL  =  M_matrLAux;
        M_rhsFull = ZeroVector( M_rhsFull.size() );

        if (M_newDomain)
            bcHdp.bdUpdate( this->mesh(), this->_feBd_u, this->_dof_u );

        bcManage( M_matrL, M_rhsFull, this->mesh(), this->_dof_u,  bcHdp , this->_feBd_u, 1.0, time );
        chrono.stop();
        std::cout << "done in " << chrono.diff() << " s." << std::endl;

        // Zero as initial guess
        M_dp = ZeroVector(  M_dp.size() );

        std::cout << "  o-  Solving system... "
                  << std::flush;
        chrono.start();

        M_linearSolverReduced.setMatrix( M_matrL );

        M_linearSolverReduced.setRecursionLevel( 1 );

        if (M_newDomain)
        {
            M_linearSolverReduced.solve(  M_dp, M_rhsFull, SolverAztec::SAME_NONZERO_PATTERN );
            M_newDomain = 0;
        }
        else
            M_linearSolverReduced.solve( M_dp, M_rhsFull, SolverAztec::SAME_PRECONDITIONER );

        chrono.stop();
        std::cout << "done in " << chrono.diff() << " s." << std::endl;

        if ( !M_linearSolverReduced.converged() )
        {
            std::cerr << "        WARNING: Solver failed to converge."
                      << std::endl;
        }
        std::cout << "        number of iterations                        = "
                  << M_linearSolverReduced.iterations() << std::endl;


        // compute pressure stress at the interface
        if ( ! bcHdz.bdUpdateDone() )
            bcHdz.bdUpdate(this->mesh(),this->_feBd_u,this->_dof_u);
        M_res = ZeroVector( M_res.size() );
        bcManageVector(M_res, this->mesh(), this->_dof_u, bcHdz, this->_feBd_u, time, 1.0 );



    }

    template <typename Mesh>
    UInt NavierStokesAleSolver<Mesh>::solveLinearSystem(Real& condEst)
    {
        M_linearSolver.setMatrix( M_matrFull );
        M_linearSolver.setRecursionLevel( 0 );

#if 1-UMFPACK_SOLVER
        if ( M_tTotalSolve < M_nUsePC * M_tThisSolve*M_tThisSolve/M_tLastSolve ) {
            solveLinearSystemOnce( false );
        } else {
            solveLinearSystemOnce( true );
            if ( !M_linearSolver.converged() ) {
                solveLinearSystemOnce( false );
            }
        }
#else
        M_linearSolver.solve( M_sol, M_rhsFull );
#endif

#if 1-UMFPACK_SOLVER
        if ( !M_linearSolver.converged() )
        {
            std::cerr << "        WARNING: Solver failed to converge."
                      << std::endl;
        }
#endif

#if PETSC_SOLVER
        condEst = M_linearSolver.condEst();
#else
        condEst = -1.;
#endif


#if 1-UMFPACK_SOLVER
        return M_linearSolver.iterations();
#else
        return 0;
#endif

    } // solveLinearSystem

    template <typename Mesh>
    void NavierStokesAleSolver<Mesh>::solveLinearSystemOnce( bool reusePC )
    {
        Chrono chrono;
        if ( reusePC ) {
            chrono.start();
            M_linearSolver.solve( M_sol, M_rhsFull,
                                  SolverAztec::SAME_PRECONDITIONER );
            chrono.stop();
            if ( M_nUsePC == 1 ) {
                M_tThisSolve = chrono.diff();
                M_tLastSolve = M_tThisSolve;
            } else {
                M_tLastSolve = M_tThisSolve;
                M_tThisSolve = chrono.diff();
            }
            M_tTotalSolve += M_tThisSolve;
            ++M_nUsePC;
        } else {
            chrono.start();
            M_linearSolver.solve( M_sol, M_rhsFull,
                                  SolverAztec::SAME_NONZERO_PATTERN );
            chrono.stop();
            M_tThisSolve = chrono.diff();
            M_tLastSolve = 2 * M_tThisSolve;
            M_tTotalSolve = M_tThisSolve;
            M_nUsePC = 1;
        }
    } // solveLinearSystemOnce


} // namespace LifeV


#endif //_NAVIERSTOKESSOLVER_H_
