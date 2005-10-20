/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Miguel A. Fernandez <miguel.fernandez@inria.fr>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2003-06-09

 Copyright (C) 2004 EPFL

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
/**
   \file NavierStokesSolver.hpp
   \author M.A. Fernandez
           C. Winkelmann
   \date 01/5/2005

   \brief This file contains a Navier-Stokes solver class which implements a
   stabilized (in space) semi-implicit (in time) scheme.
*/
#ifndef _NAVIERSTOKESSOLVER_H_
#define _NAVIERSTOKESSOLVER_H_

#define AZTEC_SOLVER 1
#define PETSC_SOLVER 0
#define UMFPACK_SOLVER 0

#include <life/lifesolver/NavierStokesHandler.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifearray/boostmatrix.hpp>

#if AZTEC_SOLVER
#include <life/lifealg/SolverAztec.hpp>
#else

#include <lifeconfig.h>

#if defined( HAVE_PETSC_H )
#include <life/lifealg/SolverPETSC.hpp>
#endif /* HAVE_PETSC_H */

#if defined ( HAVE_UMFPACK_H )
#include <life/lifealg/SolverUMFPACK.hpp>
#endif /* HAVE_UMFPACK_H */

#endif /* AZTEC_SOLVER */

#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifefem/sobolevNorms.hpp>
#include <life/lifefem/geoMap.hpp>


#include <sdStabilization.hpp>
#include <ipStabilization.hpp>


namespace LifeV
{
  /*!
    \class NavierStokesSolver

    This class implements a NavierStokes solver using different
    types of stabilizations: Interior penalty (ip), stream-line diffusion (sd)
    The resulting linear system is solved through GMRES on the full
    matrix ( u | p ).

  */
  template<typename Mesh>
  class NavierStokesSolver:
    public NavierStokesHandler<Mesh>
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
    NavierStokesSolver( const GetPot& dataFile,
                          const RefFE& refFE,
                          const QuadRule& quadRule,
                          const QuadRule& boundaryQuadRule,
                          BCHandler& bcHandler );

    //! Update the right hand side for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    void timeAdvance( source_type const& source, Real const& time );

    //! Update convective term, bc treatment and solve the linearized ns system
    void iterate( const Real& time );

    void eval( Vector& fx0, Vector& gx0, Vector x0, int status );

    void initialize( const Function& x0, Real t0=0., Real dt=0. );

    /*! Initialize when initial values for the velocity and the pressure are
      read from file
      @author Martin Prosi
    */
    void initialize( const std::string & vname );

    //! linearize convective term around given (exact) velocity function
    void linearize( const Function& betaFct ) { M_betaFct = &betaFct; }

    //! removes mean of component comp of vector x
    void removeMean( Vector& x, UInt comp=1 );


    //! Lift coefficient (y-direction)
    Real liftCoeff(const EntityFlag marker) const;


  private:

    //! Block pattern of full matrix
    MSRPatt M_fullPattern;

    //! Block pattern of full matrix
    MSRPatt M_uPattern;

#if AZTEC_SOLVER
    typedef MSRMatr<double> matrix_type;
#elif PETSC_SOLVER
    typedef BoostMatrix<boost::numeric::ublas::row_major> matrix_type;
#elif UMFPACK_SOLVER
    typedef BoostMatrix<boost::numeric::ublas::column_major> matrix_type;
#endif

    //! Matrix mass
    matrix_type M_matrMass;

    //! Matrix Stokes: rho*bdfCoeff*Vmass + mu*Vstiff + linear stabilizations
    matrix_type M_matrStokes;

    //! Matrix C: CStokes + Convective_term + nonlinear stabilizations
    matrix_type M_matrFull, M_matrFullAux;

    //! Elementary matrices and vectors
    ElemMat M_elmatC; //velocity stiffness
    ElemMat M_elmatMass; //velocity mass
    ElemMat M_elmatD;
    ElemMat M_elmatDtr;
    ElemMat M_elmatP;
    ElemVec M_elvec; // Elementary right hand side

    //! Right  hand side for the velocity
    Vector M_rhsU;

    //! Right  hand side global
    Vector M_rhsFull, M_rhsFullAux;

    //! Global solution _u and _p
    Vector M_sol;

    //! Residual vector
    Vector M_res;

#if AZTEC_SOLVER
    SolverAztec M_linearSolver;
#elif PETSC_SOLVER
    SolverPETSC M_linearSolver;
#elif UMFPACK_SOLVER
    SolverUMFPACK M_linearSolver;
#endif

    Real M_time;

    bool M_steady;

    Real M_gammaBeta;
    Real M_gammaDiv;

    const Function* M_betaFct;

    bool M_divBetaUv;
    double M_diagonalize;

    Vector M_constantPressure;

    // velocity vector for linearization of convective term
    Vector M_beta;



    SDStabilization<Mesh, Dof> M_sdStab;

    IPStabilization<Mesh, Dof> M_ipStab;

  }; // class NavierStokesSolver


  template<typename Mesh> void NavierStokesSolver<Mesh>::
  eval( Vector& fx0, Vector& /* gx0 */, Vector /* x0 */, int /* status */ )
  {
    iterate( 0.0 );
    for ( UInt iDof = 0; iDof < nDimensions*this->_dim_u ; ++iDof )
      fx0[ iDof ] = this->_u[ iDof ];

  }

  //
  //                                         IMPLEMENTATION
  //
  template<typename Mesh> NavierStokesSolver<Mesh>::
  NavierStokesSolver( const GetPot& dataFile,
			const RefFE& refFE,
			const QuadRule& quadRule,
			const QuadRule& boundaryQuadRule,
			BCHandler& bcHandler ):
    NavierStokesHandler<Mesh>( dataFile, refFE, refFE, quadRule,
                               boundaryQuadRule, quadRule, boundaryQuadRule,
                               bcHandler ),
    M_fullPattern( this->_dof_u, this->_mesh, nDimensions+1, this->patternType() ),
    M_uPattern( this->_dof_u, this->_mesh, nDimensions, this->patternType() ),
    M_matrMass( M_uPattern ),
    M_matrStokes( M_fullPattern ),
    M_matrFull( M_fullPattern ),
    M_matrFullAux( M_fullPattern ),
    M_elmatC( this->_fe_u.nbNode, nDimensions, nDimensions ),
    M_elmatMass( this->_fe_u.nbNode, nDimensions, nDimensions ),
    M_elmatD( this->_fe_u.nbNode, nDimensions+1, nDimensions ),
    M_elmatDtr( this->_fe_u.nbNode, nDimensions, nDimensions+1 ),
    M_elmatP( this->_fe_u.nbNode, nDimensions+1, nDimensions+1 ),
    M_elvec( this->_fe_u.nbNode, nDimensions ),
    M_rhsU( nDimensions*this->_dim_u ),
    M_rhsFull( ( nDimensions+1 )*this->_dim_u ),
    M_rhsFullAux( ( nDimensions+1 )*this->_dim_u ),
    M_sol( ( nDimensions+1 )*this->_dim_u ),
    M_res( ( nDimensions+1 )*this->_dim_u ),
    M_betaFct( 0 ),
    M_constantPressure( ( nDimensions+1 )*this->_dim_u ),
    M_beta( nDimensions * this->_dim_u),
    M_sdStab( dataFile, this->_mesh, this->_dof_u, this->_refFE_u, this->_Qr_u,
	      this->viscosity() ),
    M_ipStab( dataFile, this->_mesh, this->_dof_u, this->_refFE_u, this->_feBd_u, this->_Qr_u,
	      this->viscosity() )
  {
    M_steady = dataFile( "fluid/miscellaneous/steady", 1 );
    M_divBetaUv = dataFile( "fluid/discretization/div_beta_u_v", 0);
    M_diagonalize = dataFile( "fluid/discretization/diagonalize", 1.);

    // check mesh for elements with all nodes on the boundary
    UInt nLocalFaces = this->_mesh.numLocalFaces();
    UInt nFixedTets = 0;
    for ( UInt iVol = 1; iVol <= this->_mesh.numVolumes(); iVol++ )
      {
        UInt nBoundaryFaces = 0;
        for ( UInt iLocalFace=1; iLocalFace<=nLocalFaces; ++iLocalFace )
	  {
            if ( this->_mesh.isBoundaryFace( this->_mesh.localFaceId( iVol, iLocalFace )))
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
        std::cerr << "WARNING: " << nFixedTets << " of " << this->_mesh.numVolumes()
                  << " tetrahedrons have all nodes on the boundary."
                  << std::endl;
      }

#if AZTEC_SOLVER
    M_linearSolver.setOptionsFromGetPot( dataFile, "fluid/aztec" );
#elif PETSC_SOLVER
    M_linearSolver.setOptionsFromGetPot( dataFile, "fluid/petsc" );

    if ( this->bcHandler().hasOnlyEssential() && !M_diagonalize )
      {
        Real constPress = 1. / sqrt( this->_dim_u );
        for( UInt i=0; i<this->_dim_u*nDimensions; ++i )
	  {
            M_constantPressure[ i ] = 0;
	  }
        for( UInt i=this->_dim_u*nDimensions; i<this->_dim_u*(1+nDimensions); ++i )
	  {
            M_constantPressure[ i ] = constPress;
	  }
        std::vector<const Vector*> nullSpace(1);
        nullSpace[ 0 ] = &M_constantPressure;
        M_linearSolver.setNullSpace(nullSpace);
      }
#endif

    std::cout << std::endl;
    std::cout << "O-  Velocity and pressure unknowns: " << this->_dim_u     << std::endl;
    std::cout << "O-  Computing mass and Stokes matrices... " << std::flush;

    Chrono chrono;
    chrono.start();

    // Matrices initialization
    //M_u.zeros();
    M_matrMass.zeros();
    M_matrStokes.zeros();

    // Number of velocity components
    UInt nbCompU = this->_u.nbcomp();

    Real bdfCoeff = this->_bdf.bdf_u().coeff_der( 0 ) / this->_dt;

    // Elementary computation and matrix assembling
    // Loop on elements
    for ( UInt iVol = 1; iVol <= this->_mesh.numVolumes(); iVol++ )
    {
        this->_fe_u.updateFirstDeriv( this->_mesh.volumeList( iVol ) );

        M_elmatC.zero();
        M_elmatMass.zero();
        M_elmatD.zero();
        M_elmatDtr.zero();

        // stiffness strain
        stiff_strain( 2.0*this->_mu, M_elmatC, this->_fe_u );
        //stiff_div( 0.5*_fe_u.diameter(), M_elmatC, _fe_u );

        // mass
        if ( !M_steady )
        {
            mass( this->_rho*bdfCoeff, M_elmatMass, this->_fe_u, 0, 0, nDimensions );
            M_elmatC.mat() += M_elmatMass.mat();
            M_elmatMass.mat() *= ( 1./bdfCoeff );
        }

        for ( UInt iComp = 0; iComp<nbCompU; iComp++ )
        {
            for ( UInt jComp = 0; jComp<nbCompU; jComp++ )
	      // stiffness
	      assemb_mat( M_matrStokes, M_elmatC, this->_fe_u, this->_dof_u, iComp, jComp );

	    if ( !M_steady )
	      assemb_mat( M_matrMass, M_elmatMass, this->_fe_u, this->_dof_u, iComp, iComp);

            // div
            grad( iComp, 1.0, M_elmatDtr, this->_fe_u, this->_fe_u, iComp, nDimensions );
            div( iComp, -1.0, M_elmatD  , this->_fe_u, this->_fe_u, nbCompU, iComp );

            // assembling
            assemb_mat( M_matrStokes, M_elmatDtr, this->_fe_u, this->_dof_u, iComp,
                        nbCompU );
            assemb_mat( M_matrStokes, M_elmatD, this->_fe_u, this->_dof_u, nbCompU,
                        iComp );

        }
    }

    M_sol = ZeroVector( M_sol.size() );

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
    // M_matrStokes.spy( "CS.m" );


  } // NavierStokesSolver() [Constructor]

  template<typename Mesh>
  void NavierStokesSolver<Mesh>::
  timeAdvance( source_type const& source, Real const& time )
  {
    M_time = time;

    std::cout << std::endl;
    std::cout << "O== Now we are at time "<< M_time << " s." << std::endl;

    // Number of velocity components
    UInt nbCompU = this->_u.nbcomp();

    std::cout << "  o-  Updating mass term on right hand side... "
              << std::flush;

    Chrono chrono;
    chrono.start();

    // Right hand side for the velocity at time
    M_rhsU = ZeroVector( M_rhsU.size() );

    // loop on volumes: assembling source term
    for ( UInt iVol = 1; iVol<= this->_mesh.numVolumes(); ++iVol )
    {
        M_elvec.zero();
        this->_fe_u.updateJacQuadPt( this->_mesh.volumeList( iVol ) );

        for ( UInt iComp = 0; iComp<nbCompU; ++iComp )
        {
            // compute local vector
            compute_vec( source, M_elvec, this->_fe_u, M_time, iComp );

            // assemble local vector into global one
            assemb_vec( M_rhsU, M_elvec, this->_fe_u, this->_dof_u, iComp );
        }
    }

    if ( !M_steady )
      M_rhsU += M_matrMass * this->_bdf.bdf_u().time_der( this->_dt );

    M_rhsFull = ZeroVector( M_rhsFull.size() );

    for ( UInt iDof = 0; iDof<nDimensions*this->_dim_u; ++iDof )
         M_rhsFull[ iDof ] = M_rhsU[ iDof ];

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    chrono.start();

    if ( M_betaFct )
      this->uInterpolate( *M_betaFct, M_beta, M_time );
    else if ( M_steady )
      M_beta = this->_rho * this->_u; // last iteration
    else
      M_beta = this->_rho * this->_bdf.bdf_u().extrap(); // bdf extrapolation


    switch ( this->stabilization() )
      {
      case IP_STABILIZATION:
	break;
      case SD_STABILIZATION:
	std::cout << "  o-  Updating SD stabilization terms (rhs)...          "
		  << std::flush;
	M_sdStab.apply(this->_dt,M_rhsFull, M_beta, source, M_time);
	break;
      default:
	ERROR_MSG("This stabilization is not yet implemented");
      }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

  } // timeAdvance()


  template<typename Mesh>
  void NavierStokesSolver<Mesh>::iterate( const Real& time )
  {
    Chrono chrono;

    M_time = time;

    std::cout << "  o-  Copying Stokes matrix...                 "
	      << std::flush;
    chrono.start();
    M_matrFull = M_matrStokes;

    chrono.stop();

    std::cout << "done in " << chrono.diff() << " s."
              << std::endl;

    std::cout << "  o-  Updating convective volume terms...      "
              << std::flush;

    // Number of velocity components
    UInt nbCompU = this->_u.nbcomp();

    chrono.start();

    // loop on volumes
    for ( UInt iVol = 1; iVol<= this->_mesh.numVolumes(); ++iVol )
      {
        this->_fe_u.updateFirstDeriv( this->_mesh.volumeList( iVol ) ); //as updateFirstDer

        M_elmatC.zero();

        UInt eleID = this->_fe_u.currentId();
        // Non linear term, Semi-implicit approach
        // M_elvec contains the velocity values in the nodes
        for ( UInt iNode = 0 ; iNode<( UInt )this->_fe_u.nbNode ; iNode++ )
	  {
            UInt  iloc = this->_fe_u.patternFirst( iNode );
            for ( UInt iComp = 0; iComp<nbCompU; ++iComp )
	      {
                UInt ig = this->_dof_u.localToGlobal( eleID, iloc+1 )-1+iComp*this->_dim_u;
                M_elvec.vec()[ iloc+iComp*this->_fe_u.nbNode ] =  M_beta[ig];
	      }
	  }

        // Stabilising term: div u^n u v
        if ( M_divBetaUv )
	  mass_divw( 0.5, M_elvec, M_elmatC, this->_fe_u, 0, 0, nbCompU );


        // loop on components
        for ( UInt iComp = 0; iComp<nbCompU; ++iComp )
	  {
	    // compute local convective term and assembling
	    for (UInt ic = 0; ic < nbCompU; ++ic)
	      grad( ic, M_elvec, M_elmatC, this->_fe_u, this->_fe_u, iComp, iComp );
	    assemb_mat( M_matrFull, M_elmatC, this->_fe_u, this->_dof_u, iComp, iComp );
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
	M_ipStab.apply( M_matrFull, M_beta );
	break;
      case SD_STABILIZATION:
	std::cout << "  o-  Updating SD stabilization terms (matrix)...          "
		  << std::flush;
	M_sdStab.apply(this->_dt, M_matrFull, M_beta );
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

    M_rhsFullAux =  M_rhsFull;
    M_matrFullAux =  M_matrFull;

    // BC manage for the velocity
    if ( !this->bcHandler().bdUpdateDone() )
      this->bcHandler().bdUpdate( this->_mesh, this->_feBd_u, this->_dof_u );
    bcManage( M_matrFull, M_rhsFull, this->_mesh, this->_dof_u, this->bcHandler(),
	      this->_feBd_u, 1.0, M_time );

    if ( this->bcHandler().hasOnlyEssential() && M_diagonalize )
      M_matrFull.diagonalize( nDimensions*this->_dim_u, M_diagonalize, M_rhsFull, 0);

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    M_linearSolver.setMatrix( M_matrFull );



    // --------------------
    // Solving linear sytem
    // --------------------

    // set initial guess
    if ( M_steady )
      {
        // use last iteration value as initial guess
        for ( UInt i = 0; i<nDimensions*this->_dim_u; ++i )
	  M_sol[ i ] = this->_u[ i ];
      }
    else
      // use bdf based extrapolation as initial guess
      M_sol = this->_bdf.bdf_u().extrap();

    if ( this->bcHandler().hasOnlyEssential() && !M_diagonalize )
      removeMean( M_sol, 4 );


    std::cout << "  o-  Solving system...                        "
              << std::flush;
    chrono.start();
#if USE_AZTEC_SOLVER
    M_linearSolver.solve( M_sol, M_rhsFull);
    //  SolverAztec::SAME_PRECONDITIONER );
#else
    M_linearSolver.solve( M_sol, M_rhsFull );
#endif
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

#if 1-UMFPACK_SOLVER
    if ( !M_linearSolver.converged() )
      std::cerr << "        WARNING: Solver failed to converge."
		<< std::endl;

#endif

#if PETSC_SOLVER
    std::cout << "        estimated condition number (preconditioned) = "
              << M_linearSolver.condEst() << std::endl;
#endif


#if 1-UMFPACK_SOLVER
    std::cout << "        number of iterations                        = "
              << M_linearSolver.iterations() << std::endl;
#endif

    if ( this->bcHandler().hasOnlyEssential() && !M_diagonalize )
      removeMean( M_sol, 4 );

    for ( UInt iDof = 0; iDof<nDimensions*this->_dim_u; ++iDof )
      this->_u[ iDof ] = M_sol[ iDof ];

    for ( UInt iDof = 0; iDof<this->_dim_u; ++iDof )
      this->_p[ iDof ] = M_sol[ iDof+nDimensions*this->_dim_u ];

    if ( !M_steady )
      this->_bdf.bdf_u().shift_right( M_sol );

    // residual conputation for lift
    M_res = M_matrFullAux * M_sol - M_rhsFullAux;

  } // iterate()

  template <typename Mesh>
  void NavierStokesSolver<Mesh>::initialize( const Function& x0,
					       Real t0, Real dt )
  {
    ID nbComp = this->_u.nbcomp(); // Number of components of the velocity
    this->_bdf.bdf_u().initialize_unk( x0, this->_mesh, this->_refFE_u, this->_fe_u, this->_dof_u, t0,
                                 dt, nbComp+1 );

    // initialize M_sol with the first element in bdf_u.unk (=last value)
    M_sol = *( this->_bdf.bdf_u().unk().begin() );
    if ( this->bcHandler().hasOnlyEssential() && !M_diagonalize )
      {
        removeMean( M_sol, 4 );
      }
    for ( UInt iDof = 0; iDof<nDimensions*this->_dim_u; ++iDof )
      {
        this->_u[ iDof ] = M_sol[ iDof ];
      }
    for ( UInt iDof = 0; iDof<this->_dim_u; ++iDof )
      {
        this->_p[ iDof ] = M_sol[ iDof+nDimensions*this->_dim_u ];
      }

    // write initial values (for debugging only)
    const std::vector<Vector>& unk = this->_bdf.bdf_u().unk();
    for ( int iUnk=unk.size()-1; iUnk>=0; --iUnk )
      {
        for ( UInt iDof=0; iDof<nDimensions*this->_dim_u; ++iDof )
	  {
            this->_u[ iDof ] = unk[ iUnk ][ iDof ];
	  }
        for ( UInt iDof = 0; iDof<this->_dim_u; ++iDof )
	  {
            this->_p[ iDof ] = unk[ iUnk ][ iDof+nDimensions*this->_dim_u ];
	  }
        //this->postProcess();
      }
  } // initialize()

  /*! Initialize when initial values for the velocity and the pressure are read
    from file
    @author Martin Prosi
  */
  template <typename Mesh>
  void NavierStokesSolver<Mesh>::initialize( const std::string & vname )
  {

    std::fstream resfile( vname.c_str(), std::ios::in | std::ios::binary );
    if ( resfile.fail() )
      {
        std::cerr << " Error in initialization: File not found or locked"
                  << std::endl;
        abort();
      }

    //  resfile.read( ( char* ) & M_sol( 0 ), M_sol.size() * sizeof( double ) );
    // resfile.close();

    for ( UInt iDof = 0; iDof<nDimensions*this->_dim_u; ++iDof )
      {
        resfile >> M_sol[ iDof ]  ;
        this->_u[ iDof ] = M_sol[ iDof ];
      }
    for ( UInt iDof = 0; iDof<this->_dim_u; ++iDof )
      {
	resfile >> M_sol[iDof+nDimensions*this->_dim_u  ]  ;
        this->_p[ iDof ] = M_sol[ iDof+nDimensions*this->_dim_u ];
      }

    this->_bdf.bdf_u().initialize_unk( M_sol );
    //_bdf.bdf_p().initialize_unk( _p );

    this->_bdf.bdf_u().showMe();
    //_bdf.bdf_p().showMe();

  }



  template <typename Mesh>
  void NavierStokesSolver<Mesh>::removeMean( Vector& x, UInt comp )
  {
    Real sum1 = 0.;
    Real sum0 = 0.;
    for ( UInt iVol = 1; iVol <= this->_mesh.numVolumes(); iVol++ )
      {
        this->_fe_p.updateFirstDeriv( this->_mesh.volumeList( iVol ) );
        sum1 += elem_integral( x, this->_fe_p, this->_dof_p, comp );
        sum0 += this->_fe_p.measure();
      }
    Real mean = sum1/sum0;
    for ( UInt iNode = (comp-1)*this->_dim_u; iNode<comp*this->_dim_u; ++iNode )
      {
        x[ iNode ] -= mean;
      }

  } // removeMean()



  template <typename Mesh>
  Real NavierStokesSolver<Mesh>::liftCoeff(const EntityFlag marker) const
  {

    Real liftCoeff = 0;

    for (ID i=1; i <= this->_dim_u; ++i)
      if (this->_mesh.pointList( i ).marker() == marker)
	liftCoeff += M_res[ this->_dim_u + i-1 ];

    return liftCoeff;

  }







} // namespace LifeV

#endif //_NAVIERSTOKESSOLVER_H_
