/*
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
/*!
  \file NavierStokesAleSolverCT.h
  \author M.A. Fernandez
  \date 02/2005
  \version 1.0

  \brief This file contains a Navier-Stokes solver allowing moving domains (via an
  ALE formulation). This class implements a time Chorin-Teman semi-implicit scheme with stabilized
  finite elements in space (SUPG). Tested with P1/P1 and Q1/Q1.

*/

#ifndef _NAVIERSTOKESALESOLVERCT_HH
#define _NAVIERSTOKESALESOLVERCT_HH

#include <life/lifesolver/NavierStokesHandler.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifefem/assembGeneric.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifecore/chrono.hpp>
//#include <life/lifefem/meshMotion.hpp>
//#include <life/lifealg/SolverAztec.hpp>
#include <life/lifesolver/sdStabilization.hpp>


namespace LifeV
{

  /*!
    \class NavierStokesAleSolverCT


  */
  template <typename Mesh>
  class NavierStokesAleSolverCT:
    public NavierStokesHandler<Mesh>
  {

  public:

    typedef typename NavierStokesHandler<Mesh>::Function Function;
    typedef typename NavierStokesHandler<Mesh>::source_type source_type;

    //! Constructor
    /*!
      \param dataFile GetPot data file
      \param refFE_u reference FE for the velocity
      \param Qr_u volumic quadrature rule for the velocity
      \param bdQr_u surface quadrature rule for the velocity
      \param bcHu boundary conditions for the velocity
      \param bcHp boundary conditions for the pressure
      \param BCh_mesh boundary conditions for the motion harmonic extension
    */
    NavierStokesAleSolverCT( const GetPot& dataFile,
                             const RefFE& refFE_u,
			     const QuadRule& Qr_u,
                             const QuadRule& bdQr_u,
			     BCHandler& bcHu,
			     BCHandler& bcHp,
                             BCHandler& bcExtOper );


    void timeAdvance(source_type const& source, Real const& time );
    void iterate(const Real& time);
    void iterateLin(const Real time, BCHandler& bcHdu, BCHandler& bcHdp, const bool shapeTerms=1);
    void iteratePressureLin(const Real time, BCHandler& bcHdp);
    void updateMesh( const Real& time );
    void dwUpdate( const Real& time );

    void iterateVelocity(const Real& time);
    void iteratePressure(const Real& time);

    Vector& w() {return M_w;}
    Vector& dw() {return M_dw;}
    const Vector& d() {return M_d;}
    Vector& dp() {return M_dp;}
    void updateForce(BCHandler& bcHd);
    void updateLinForce(BCHandler& bcHdz);
    Vector& residual() {return M_force;}

  private:


    //! Block pattern of M_u
    MSRPatt M_pattM_u_block;

    //! Pattern for M
    MixedPattern<3, 3, MSRPatt> M_pattM_u;

    //! pattern velocity and pressure matrix
    MSRPatt M_pattC_u, M_pattC_p;

    MixedMatr<3, 3, MSRPatt, double> M_M_u;

    MSRMatr<double> M_C_u, M_C_p, M_C_uWithOutBC, M_C_pWithOutBC;

    //! Elementary matrices and vectors
    ElemMat M_elmatM_u, M_elmatC_u, M_elmatC_p;
    ElemVec M_elvec_u, M_elvec_p, M_uLoc, M_betaLoc, M_pLoc, M_wLoc;

    //! Right  hand  side for the velocity
    Vector M_f_u, M_f_p, M_du, M_dp, M_dw, M_residual_u,
      M_residual_du,M_f_uWithOutBC, M_f_pWithOutBC, M_force, M_f_uAux;

    //! linear solvers: velocity and pressure
    SolverAztec M_linearSolver_u, M_linearSolver_p;


    HarmonicExtension M_extOper;
    const Vector& M_d;

    //! Backward velocity and mesh velocity
    Vector M_un, M_w, M_beta, M_dOld;

    BCHandler& M_bcH_p;

    bool M_firstTimeStep;
    bool M_newDomain;

    SDStabilization<Mesh, Dof> M_sdStab;

    void M_computeVelocityMassMatrix();

    void M_solveVelocitySystem( const Real& time );
    void M_computePressureRHS(const Vector& tilde_u);
    void M_computePressureMatrix();
    void M_solvePressureSystem( const Real& time );

    UInt M_nUsePC_u;
    Real M_tTotalSolve_u;
    Real M_tLastSolve_u;
    Real M_tThisSolve_u;

    UInt M_nUsePC_p;
    Real M_tTotalSolve_p;
    Real M_tLastSolve_p;
    Real M_tThisSolve_p;

    /*! solve linear system
      @return number of iterations, zero if not applicable
      @param condEst condition estimate is returned here, -1 if not applicable
    */
    UInt solveLinearSystem_u();
    void solveLinearSystemOnce_u( bool reusePC );
    UInt solveLinearSystem_p();
    void solveLinearSystemOnce_p( bool reusePC );

    Real M_time;


  };


  //
  //                                         IMPLEMENTATION
  //
  template <typename Mesh>
  NavierStokesAleSolverCT<Mesh>::
  NavierStokesAleSolverCT( const GetPot& dataFile, const RefFE& refFE_u, const QuadRule& Qr_u,
			   const QuadRule& bdQr_u, BCHandler& bcHu,
			   BCHandler& bcHp, BCHandler& bcExtOper ) :
    NavierStokesHandler<Mesh>( dataFile, refFE_u, refFE_u, Qr_u, bdQr_u, Qr_u, bdQr_u, bcHu ),
    M_pattM_u_block( this->_dof_u ),
    M_pattM_u( M_pattM_u_block, "diag" ),
    M_pattC_u( this->_dof_u, this->_fe_u.nbCoor ),
    M_pattC_p( this->_dof_u, 1 ),
    M_M_u( M_pattM_u ),
    M_C_u( M_pattC_u ),
    M_C_p( M_pattC_p ),
    M_C_uWithOutBC( M_pattC_u ),
    M_C_pWithOutBC( M_pattC_p ),
    M_elmatM_u( this->_fe_u.nbNode, this->_fe_u.nbCoor, this->_fe_u.nbCoor),
    M_elmatC_u( this->_fe_u.nbNode, this->_fe_u.nbCoor, this->_fe_u.nbCoor),
    M_elmatC_p( this->_fe_u.nbNode, 1, 1),
    M_elvec_u(  this->_fe_u.nbNode, this->_fe_u.nbCoor ),
    M_elvec_p(  this->_fe_u.nbNode, 1),
    M_uLoc( this->_fe_u.nbNode, this->_fe_u.nbCoor ),
    M_betaLoc( this->_fe_u.nbNode, this->_fe_u.nbCoor ),
    M_pLoc( this->_fe_u.nbNode, 1 ),
    M_wLoc( this->_fe_u.nbNode,  this->_fe_u.nbCoor ),
    M_f_u( this->_fe_u.nbCoor*this->_dim_u ),
    M_f_p( this->_dim_p ),
    M_du( this->_fe_u.nbCoor*this->_dim_u ),
    M_dp( this->_dim_p ),
    M_dw(  this->_fe_u.nbCoor*this->_dim_u ),
    M_residual_u(  this->_fe_u.nbCoor*this->_dim_u ),
    M_residual_du(  this->_fe_u.nbCoor*this->_dim_u ),
    M_f_uWithOutBC( this->_fe_u.nbCoor*this->_dim_u ),
    M_f_pWithOutBC( this->_dim_p ),
    M_force(this->_fe_u.nbCoor*this->_dim_u ),
    M_f_uAux(this->_fe_u.nbCoor*this->_dim_u ),
    M_extOper(  dataFile,
		this->_mesh, 1.0,
		Qr_u,
		bdQr_u,
		bcExtOper ),
    M_d(  M_extOper.getDisplacement()  ),
    M_un(  this->_fe_u.nbCoor*this->_dim_u ),
    M_w(  this->_fe_u.nbCoor*this->_dim_u ),
    M_beta( this->_fe_u.nbCoor*this->_dim_u ),
    M_dOld( this->_fe_u.nbCoor*this->_dim_u ),
    M_bcH_p(bcHp),
    M_firstTimeStep(1),
    M_newDomain(0),
    M_sdStab( dataFile, this->_mesh, this->_dof_u, this->_refFE_u, this->_Qr_u, this->viscosity() ),
    M_nUsePC_u( 1 ),
    M_tTotalSolve_u( 0 ),
    M_tLastSolve_u( 1 ),
    M_tThisSolve_u( 1 ),
    M_nUsePC_p( 1 ),
    M_tTotalSolve_p( 0 ),
    M_tLastSolve_p( 1 ),
    M_tThisSolve_p( 1 ),
    M_time(0)
  {


    std::cout << std::endl;
    std::cout << "O-  Pressure unknowns: " << this->_dim_p << std::endl;
    std::cout << "O-  Velocity unknowns: " << this->_dim_u << std::endl << std::endl;


    M_linearSolver_u.setOptionsFromGetPot( dataFile, "fluid/aztec_u" );


    M_linearSolver_p.setOptionsFromGetPot( dataFile, "fluid/aztec_p" );

    // initialize pressure
    this->_p = ZeroVector( this->_p.size() );
    M_residual_u = ZeroVector( M_residual_u.size() );
    M_residual_du = ZeroVector( M_residual_du.size() );
    M_force = ZeroVector( M_force.size() );

  }


  template <typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::updateMesh( const Real& time )
  {


    Chrono chrono;

    std::cout << "\n  o-  Updating fluid mesh... ";

    // Updating mesh displacement and velocity
    this->M_extOper.updateExtension( this->_mesh, time );

    // Updating mesh points
    this->_mesh.moveMesh( M_d );

    Real dti = 1.0 / this->_dt;
    M_w = ( M_d - M_dOld ) * dti;

    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    std::cout << "   norm( df^{n+1} )  = " << norm_inf( M_d ) << std::endl;
    std::cout << "   norm( wf^{n+1} )  = " << norm_inf( M_w ) << std::endl;

  }

  template <typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::dwUpdate( const Real& time )
  {


    Chrono chrono;

    std::cout << "\n  o-  Updating linearized mesh velocity... ";

    // Updating mesh displacement and velocity
    this->M_extOper.updateExtension( this->_mesh, time );

    Real dti = 1.0 / this->_dt;
    M_dw = M_d * dti;
    std::cout << "   norm( M_dw )  = " << norm_inf( M_dw ) << std::endl;
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

  }



  template <typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::timeAdvance(source_type const& source, Real const& time )
  {


    M_time = time;

    std::cout << std::endl;
    std::cout << "O== FLUID: Now we are at time " << time << " s." << std::endl;

    UInt ig;

    Chrono chrono;
    chrono.start();

    // compute mass matrix at the first time step
    if ( M_firstTimeStep )
      {
	std::cout << "  F-  Computing mass matrix and updating mass and source term on right hand side... "
		  << std::flush;
	this->M_computeVelocityMassMatrix();
	M_firstTimeStep = 0;
      }
    else
      std::cout << "  F-  Updating mass and source term on right hand side... " << std::flush;

    M_f_u = ZeroVector( M_f_u.size() );

    // Loop on elements
    for ( UInt i = 1; i <= this->_mesh.numVolumes(); i++ )
      {

	this->_fe_u.updateFirstDeriv( this->_mesh.volumeList( i ) );

	M_elvec_u.zero();

	// updating pressure  \tilde p^{n} on the current element
	for ( UInt k = 0; k < ( UInt ) this->_fe_u.nbNode; ++k )
	  {
	    ig = this->_dof_u.localToGlobal( i, k + 1 ) - 1;
	    M_pLoc[ k ] = this->_p[ig];
	  }

	// we add the term - ( \grad p^{n}, v)
	for (UInt ic=0; ic< (UInt)this->_fe_u.nbCoor; ++ic)
	  {
	    source_gradpv(-1.0, M_pLoc, M_elvec_u, this->_fe_u, this->_fe_u, ic);

	    // assembling into global vector
	    assemb_vec(M_f_u, M_elvec_u, this->_fe_u, this->_dof_u, ic);
	  }
      }

    // finaly we add the mass term: (rho / dt) * u^n
    M_f_u += M_M_u * this->_u;
    M_f_uWithOutBC= M_f_u;

    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    // updating old displacement
    M_dOld = M_d;
    M_un = this->_u;
  }

  // implicit: solves velocity and pressure system
  template <typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::iterate( const Real& time )
  {

    M_solveVelocitySystem( time );

    Chrono chrono;
    std::cout << "  o-  Updating pressure right hand side... ";
    chrono.start();
    M_computePressureRHS( this->_u );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    M_solvePressureSystem( time );
  }


  // semi-implicit: solves velocity system and updates pressure rhs and matrix (without bc)
  template <typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::iterateVelocity( const Real& time )
  {

    // solve velocity system
    M_solveVelocitySystem( time );

    // compute pressure system: rhs and matrix (without bc)
    Chrono chrono;
    std::cout << "  o-  Updating pressure right hand side... ";
    chrono.start();
    M_computePressureRHS( this->_u );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    // rhs constant on pressure inner iterations
    M_f_pWithOutBC= M_f_p;

    M_computePressureMatrix();

    // new pressure matrix
    M_newDomain = 1;



  }

  // semi-implicit: solves pressure system
  template <typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::iteratePressure( const Real& time )
  {

    M_f_p = M_f_pWithOutBC; // computed from iterate velocity once per time step

    Chrono chrono;

    // recover matrix without
    M_C_p = M_C_pWithOutBC;

    std::cout << "  o-  Applying pressure boundary conditions... ";
    chrono.start();
    if ( ! M_bcH_p.bdUpdateDone() )
      M_bcH_p.bdUpdate(this->_mesh,this->_feBd_u,this->_dof_u);
    bcManage(M_C_p,M_f_p,this->_mesh,this->_dof_u,M_bcH_p,this->_feBd_u,1.0,time);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    std::cout << "  o-  Solving pressure system... "<< std::flush;
    chrono.start();
    if ( M_newDomain )
      {

	this->_p = ZeroVector( this->_p.size() );
	UInt iter = solveLinearSystem_p();
	M_newDomain = 0;
      }
    else
      M_linearSolver_p.solve( this->_p, M_f_p, SolverAztec::SAME_PRECONDITIONER);

    chrono.stop();

    std::cout << "done in " << chrono.diff() << " s." << std::endl;



  }

  // linearized version: for FSI, implicit scheme solved via newton
  template <typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::iterateLin(const Real time, BCHandler& bcHdu, BCHandler& bcHdp, const bool shapeTerms)
  {


    Chrono chrono;
    UInt  nbCompU =  this->_fe_u.nbCoor;
    UInt iloc,ig,icloc;

    ElemVec dLoc( this->_fe_u.nbNode,  this->_fe_u.nbCoor );
    ElemVec dwLoc( this->_fe_u.nbNode,  this->_fe_u.nbCoor );
    ElemVec unLoc( this->_fe_u.nbNode,  this->_fe_u.nbCoor );

    //
    // linearized velocity system
    //

    M_f_uAux = ZeroVector( M_f_uAux.size() );



    if (shapeTerms) {
      /*

      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      No taking into account source terms
      i.e f = 0. To be modified
      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      */

      std::cout << "    F-LINEAR-  Updating right hand side (shape derivative terms)... ";


      chrono.start();

      // Elementary computation and right hand-side assembling
      // Loop on elements
      for ( UInt i = 1; i <= this->_mesh.numVolumes(); ++i )
	{
	  this->_fe_u.updateFirstDeriv( this->_mesh.volumeList( i ) );


	  M_elvec_u.zero();

	  for ( UInt k = 0 ; k < ( UInt ) this->_fe_u.nbNode ; ++k )
	    {
	      iloc = this->_fe_u.patternFirst( k );
	      for ( UInt ic = 0; ic < nbCompU; ++ic )
		{
		  ig = this->_dof_u.localToGlobal( i, iloc + 1 ) - 1 + ic * this->_dim_u;
		  icloc = iloc + ic * this->_fe_u.nbNode ;
		  M_betaLoc.vec()[ icloc  ] = M_un( ig ) - M_w( ig );  // u^n - w^k local
		  M_wLoc.vec()   [ icloc  ] = M_w( ig );               // w^k local
		  M_uLoc.vec()   [ icloc  ] = this->_u( ig );          // u^k local
		  dLoc.vec()   [ icloc  ] = M_d( ig );               // d local
		  dwLoc.vec()  [ icloc  ] = M_dw( ig );              // dw local
		  unLoc.vec()  [ icloc  ] = M_un( ig );              // u^n local
		}
	      ig = this->_dof_u.localToGlobal( i, iloc + 1 ) - 1;
	      M_pLoc[ iloc ] = 0.0;  // p^k local
	    }

	  //  - \rho ( -\grad w^k :[I\div d - (\grad d)^T] u^k + ( u^n-w^k )^T[I\div d - (\grad d)^T] (\grad u^k)^T , v  )
	  source_mass1( -this->_rho, M_uLoc, M_wLoc, M_betaLoc, dLoc, M_elvec_u, this->_fe_u );

	  //  + \rho * ( \grad u^k dw, v  )
	  source_mass2( this->_rho, M_uLoc, dwLoc, M_elvec_u, this->_fe_u );

	  //  - \rho/2 ( \nabla u^n:[2 * I\div d - (\grad d)^T]  u^k , v  )
	  source_mass3( -0.5*this->_rho, unLoc, M_uLoc, dLoc, M_elvec_u, this->_fe_u );

	  //  - ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  )
	  source_stress( -1.0, this->_mu, M_uLoc, M_pLoc, dLoc, M_elvec_u, this->_fe_u, this->_fe_u );

	  // + \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )
	  source_stress2( this->_mu, M_uLoc, dLoc, M_elvec_u, this->_fe_u );

	  // loop on velocity components
	  for ( UInt ic = 0; ic < nbCompU ; ic++ )
	    {
	      // assembling velocity right hand side
	      assemb_vec( M_f_uAux, M_elvec_u, this->_fe_u, this->_dof_u, ic );
	    }
	}

    }



    // matrix known
    M_C_u = M_C_uWithOutBC;


    M_f_u = M_f_uAux;
    std::cout << "   norm( M_f_du )  = " << norm_inf( M_f_uAux ) << std::endl;


    std::cout << "  o-  Applying linearized velocity boundary conditions... ";
    chrono.start();
    if ( ! bcHdu.bdUpdateDone() )
      bcHdu.bdUpdate(this->_mesh,this->_feBd_u,this->_dof_u);
    bcManage(M_C_u, M_f_u,this->_mesh, this->_dof_u, bcHdu,this->_feBd_u, 1.0, time);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    M_du = ZeroVector( M_du.size() );

    std::cout << "  o-  Solving linearized velocity system... "<< std::flush;
    chrono.start();
    M_linearSolver_u.solve( M_du, M_f_u, SolverAztec::SAME_PRECONDITIONER);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    M_residual_du =   M_f_uAux - M_C_uWithOutBC * M_du;


    //
    // linearized pressure system
    //
    std::cout << "  o-  Updating linearized pressure right hand side... ";
    chrono.start();
    M_computePressureRHS( M_du );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;



    if (shapeTerms)
      {
	/*

	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	No taking into account source terms
	i.e f = 0. To be modified
	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	*/

	std::cout << "    F-LINEAR-  Updating pressure rhs (shape derivative terms)... ";
	UInt iloc,ig,icloc;

	chrono.start();

	// Elementary computation and right hand-side assembling
	// Loop on elements
	for ( UInt i = 1; i <= this->_mesh.numVolumes(); ++i )
	  {
	    this->_fe_u.updateFirstDeriv( this->_mesh.volumeList( i ) );


	    M_elvec_p.zero();

	    for ( UInt k = 0 ; k < ( UInt ) this->_fe_u.nbNode ; ++k )
	      {
		iloc = this->_fe_u.patternFirst( k );
		for ( UInt ic = 0; ic < nbCompU; ++ic )
		  {
		    ig = this->_dof_u.localToGlobal( i, iloc + 1 ) - 1 + ic * this->_dim_u;
		    icloc = iloc + ic * this->_fe_u.nbNode ;
		    dLoc.vec()   [ icloc  ] = M_d( ig );
		  }
		ig = this->_dof_u.localToGlobal( i, iloc + 1 ) - 1;
		M_pLoc[ iloc ] = this->_p(ig);
	      }

	    //  - (  [I\div d - (\grad d)^T - \grad d]\grap p^k, \grad q  )
	    source_press2( -1.0, M_pLoc, dLoc, M_elvec_p, this->_fe_u);
	    assemb_vec( M_f_p, M_elvec_p, this->_fe_u, this->_dof_u);

	  }

      }



    // matrix known
    M_C_p = M_C_pWithOutBC;

    std::cout << "  o-  Applying linearized pressure boundary conditions... ";
    chrono.start();
    if ( ! bcHdp.bdUpdateDone() )
      bcHdp.bdUpdate(this->_mesh,this->_feBd_u,this->_dof_u);
    bcManage(M_C_p, M_f_p, this->_mesh, this->_dof_u, bcHdp, this->_feBd_u, 1.0, time);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    M_linearSolver_p.setRecursionLevel( 0 );
    M_dp = ZeroVector( M_dp.size() );

    std::cout << "  o-  Solving linearized pressure system... "<< std::flush;
    chrono.start();
    M_linearSolver_p.solve( M_dp, M_f_p, SolverAztec::SAME_PRECONDITIONER);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

  }


  // linearized reduced version: for FSI, implicit scheme solved via quasi-newton
  template <typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::iteratePressureLin(const Real time, BCHandler& bcHdp)
  {


    Chrono chrono;

    // no viscous stress variations
    M_residual_du =  ZeroVector(  M_residual_du.size()  );

    // zero rhs
    M_f_p =  ZeroVector(  M_f_p.size()  );
    M_C_p = M_C_pWithOutBC;

    std::cout << "  o-  Applying linearized pressure boundary conditions... ";
    chrono.start();
    if ( ! bcHdp.bdUpdateDone() )
      bcHdp.bdUpdate(this->_mesh,this->_feBd_u,this->_dof_u);
    bcManage(M_C_p, M_f_p, this->_mesh, this->_dof_u, bcHdp, this->_feBd_u, 1.0, time);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    M_linearSolver_p.setRecursionLevel( 0 );
    M_dp = ZeroVector( M_dp.size() );

    std::cout << "  o-  Solving linearized pressure system... "<< std::flush;
    chrono.start();
    M_linearSolver_p.solve( M_dp, M_f_p, SolverAztec::SAME_PRECONDITIONER);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

  }


  template<typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::updateForce(BCHandler& bcHd)
  {

    if ( ! bcHd.bdUpdateDone() )
      bcHd.bdUpdate(this->_mesh, this->_feBd_u, this->_dof_u);

    // viscous stress
    M_force = M_residual_u;

    std::cout << "   norm( M_force )  = " << norm_inf( M_force ) << std::endl;
    // pressure stress
    bcManageVector(M_force, this->_mesh, this->_dof_u, bcHd, this->_feBd_u, 0.0, 1.0 );
    std::cout << "   norm( M_force )  = " << norm_inf( M_force ) << std::endl;


    std::ofstream forcef("force.bb");
    forcef << "%  time = " << M_time << std::endl;
    forcef << " 3 1 ";
    forcef << this->_dim_u;
    forcef << " 2\n";
    Real s=0;
    for (UInt i=0; i < this->_dim_u; ++i)
      {
	for (UInt j=0; j < 3; ++j)
	  s+=  M_force[ j*this->_dim_u+i]*M_force[ j*this->_dim_u+i];
	forcef << sqrt( s ) << "\n";
	s=0;
      }
    forcef.close();


  }

  template<typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::updateLinForce(BCHandler& bcHdz)
  {

    if ( ! bcHdz.bdUpdateDone() )
      bcHdz.bdUpdate(this->_mesh, this->_feBd_u, this->_dof_u);

    M_force = M_residual_du;


    // pressure stress
    bcManageVector(M_force, this->_mesh, this->_dof_u, bcHdz, this->_feBd_u, 0.0, 1.0 );
  }

  // Compute the velocity mas Matrix at the initialization
  template <typename Mesh> void NavierStokesAleSolverCT<Mesh>::M_computeVelocityMassMatrix()
  {

    // Initializing mass matrix
    M_M_u.zeros();

    Real dti = 1.0 / this->_dt;

    // loop on volumes: assembling mass term
    for ( UInt i = 1; i <= this->_mesh.numVolumes(); ++i )
      {
        this->_fe_u.updateJac( this->_mesh.volumeList( i ) );
        M_elmatM_u.zero();
        mass( this->_rho * dti, M_elmatM_u, this->_fe_u, 0, 0,  this->_fe_u.nbCoor );

	for ( UInt ic = 0; ic < (UInt)this->_fe_u.nbCoor; ++ic )
	  assemb_mat( M_M_u, M_elmatM_u, this->_fe_u, this->_dof_u, ic, ic );
      }

  }


  // Compute and solves velocity system
  template <typename Mesh> void NavierStokesAleSolverCT<Mesh>::M_solveVelocitySystem( const Real& time)
  {


    Chrono chrono;

    std::cout << "  o-  Updating velocity matrix... ";

    M_beta = this->_rho * ( M_un - M_w );

    chrono.start();

    //initialize matrices
    M_M_u.zeros();
    M_C_u.zeros();				\

    Real dti = 1.0 / this->_dt;
    UInt ig;

    // Loop on elements
    for ( UInt i = 1; i <= this->_mesh.numVolumes(); i++ )
      {

	// assem_mat_mixed
	this->_fe_u.updateFirstDerivQuadPt( this->_mesh.volumeList( i ) );

	// initialization of elementary matrices
	M_elmatM_u.zero();
	M_elmatC_u.zero();

	// mass
	mass( this->_rho * dti, M_elmatM_u, this->_fe_u, 0, 0,  this->_fe_u.nbCoor );

	// stiffness strain
	stiff_strain( 2.0 * this->_mu, M_elmatC_u, this->_fe_u );
	M_elmatC_u.mat() += M_elmatM_u.mat();

	// Non linear term, Semi-implicit approach
	// (\tilde u^{n} - w^{n+1})  on the current element
	for ( UInt k = 0; k < ( UInt ) this->_fe_u.nbNode; ++k )
	  {
	    for ( UInt ic = 0; ic < (UInt)this->_fe_u.nbCoor; ++ic )
	      {
		ig = this->_dof_u.localToGlobal( i, k + 1 ) - 1 + ic * this->_dim_u;
		M_betaLoc[ k + ic * this->_fe_u.nbNode ] = M_beta[ig];
		M_uLoc[ k + ic * this->_fe_u.nbNode ] = M_un[ig];
		M_wLoc[ k + ic * this->_fe_u.nbNode ] = M_w[ig];
	      }
	  }

	// ALE term: -  \rho * (div w u, v)
	mass_divw( -this->_rho, M_wLoc, M_elmatC_u, this->_fe_u, 0, 0, this->_fe_u.nbCoor );

	// ALE term: 0.5 * \rho * (div u^n u, v)
	mass_divw( 0.5 * this->_rho, M_uLoc, M_elmatC_u, this->_fe_u, 0, 0, this->_fe_u.nbCoor );

        // convective term (\tilde u^{n} - w^{n+1})
        for ( UInt ic = 0;ic < (UInt)this->_fe_u.nbCoor;ic++ )
	    for ( UInt jc = 0;jc < (UInt)this->_fe_u.nbCoor;jc++ )
		grad( jc, M_betaLoc, M_elmatC_u, this->_fe_u, this->_fe_u, ic, ic );


	// assembling
	for ( UInt ic = 0;ic < (UInt)this->_fe_u.nbCoor;ic++ )
	  {
	    for ( UInt jc = 0;jc < (UInt)this->_fe_u.nbCoor;jc++ )
		assemb_mat( M_C_u, M_elmatC_u, this->_fe_u, this->_dof_u, ic, jc );

	    // mass
	    assemb_mat( M_M_u, M_elmatM_u, this->_fe_u, this->_dof_u, ic, ic );
	  }
      }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;


    chrono.start();

    switch ( this->stabilization() )
      {
      case SD_STABILIZATION:
	std::cout << "  o-  Updating SD stabilization terms (matrix)...          "
		  << std::flush;
	M_sdStab.applyCT(this->_dt, M_C_u, M_beta );

	// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	// source has to be puted here
	// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	//
	//	std::cout << "  o-  Updating SD stabilization terms (rhs)...          "
	//	  << std::flush;
	//M_sdStab.applyCT(this->_dt,M_rhsFullAux, M_beta, source, M_time);
	break;
      default:
	ERROR_MSG("This stabilization is not yet implemented");
      }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    // used for residual and  linearized iterations
    M_C_uWithOutBC = M_C_u;
    M_f_u = M_f_uWithOutBC;

    std::cout << "  o-  Applying velocity boundary conditions... ";
    chrono.start();
    this->bcHandler().bdUpdate(this->_mesh,this->_feBd_u,this->_dof_u);
    bcManage(M_C_u, M_f_u,this->_mesh, this->_dof_u, this->bcHandler(),this->_feBd_u, 1.0, time);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;


    this->_u = ZeroVector( this->_u.size() );

    std::cout << "  o-  Solving velocity system... "<< std::flush;
    chrono.start();
    UInt iter = solveLinearSystem_u();
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    if ( iter > 0 )
      std::cout << "        number of iterations                        = "
		<< iter << std::endl;

    // residual: viscous stresses
    M_residual_u = M_f_uWithOutBC - M_C_uWithOutBC * this->_u;


  }



  // Computes pressure system rhs from a given velocity field (nodal)
  template <typename Mesh> void NavierStokesAleSolverCT<Mesh>::M_computePressureRHS(const Vector& tilde_u)
  {


    Real dti = 1.0 / this->_dt;
    UInt ig;

    M_f_p = ZeroVector( M_f_p.size() );

    // Loop on elements
    for ( UInt i = 1; i <= this->_mesh.numVolumes(); i++ )
      {
	this->_fe_u.updateFirstDeriv( this->_mesh.volumeList( i ) );
	M_elvec_p.zero();

	// updating velocity  \tilde u^{n+1} on the current element
	for ( UInt k = 0; k < ( UInt ) this->_fe_u.nbNode; ++k )
	  {
	    for ( UInt ic = 0; ic <  (UInt)this->_fe_u.nbCoor; ++ic )
	      {
		ig = this->_dof_u.localToGlobal( i, k + 1 ) - 1 + ic * this->_dim_u;
		M_uLoc[ k + ic * this->_fe_u.nbNode ] = tilde_u[ig];
	      }
	  }

	// we add the term - (1/dt) * ( \div \tilde u^{n+1}, q)
	source_divuq( -1.0*dti, M_uLoc, M_elvec_p, this->_fe_u, this->_fe_u);
	assemb_vec(M_f_p, M_elvec_p, this->_fe_u, this->_dof_u);
      }

  }


  // Compute the pressure matrix
  template <typename Mesh> void NavierStokesAleSolverCT<Mesh>::M_computePressureMatrix( )
  {

    Chrono chrono;

    std::cout << "  o-  Updating pressure matrix... ";
    chrono.start();

    M_C_p.zeros();

    // Loop on elements
    for ( UInt i = 1; i <= this->_mesh.numVolumes(); i++ )
      {
	this->_fe_u.updateFirstDeriv( this->_mesh.volumeList( i ) );

	M_elmatC_p.zero();

	stiff(1.0, M_elmatC_p, this->_fe_u );

	assemb_mat( M_C_p, M_elmatC_p, this->_fe_u, this->_dof_u, 0, 0);

      }
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    // used for linearized iterations (Newton)
    M_C_pWithOutBC= M_C_p;

  }




  // Compute the pressure matrix and solves the system
  template <typename Mesh> void NavierStokesAleSolverCT<Mesh>::M_solvePressureSystem( const Real& time )
  {

    M_computePressureMatrix();

    Chrono chrono;

    std::cout << "  o-  Applying pressure boundary conditions... ";
    chrono.start();
    if ( ! M_bcH_p.bdUpdateDone() )
      M_bcH_p.bdUpdate(this->_mesh,this->_feBd_u,this->_dof_u);
    bcManage(M_C_p,M_f_p,this->_mesh,this->_dof_u,M_bcH_p,this->_feBd_u,1.0,time);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;


    std::cout << "  o-  Solving pressure system... "<< std::flush;
    chrono.start();
    this->_p = ZeroVector( this->_p.size() );
    UInt iter = solveLinearSystem_p();
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    if ( iter > 0 )
      std::cout << "        number of iterations                        = "
		<< iter << std::endl;

  }


  template <typename Mesh>
  UInt NavierStokesAleSolverCT<Mesh>::solveLinearSystem_u()
  {

    M_linearSolver_u.setMatrix( M_C_u );
    M_linearSolver_u.setRecursionLevel( 0 );

    if ( M_tTotalSolve_u < M_nUsePC_u * M_tThisSolve_u*M_tThisSolve_u/M_tLastSolve_u )
      solveLinearSystemOnce_u( false );
    else
      {
        solveLinearSystemOnce_u( true );
        if ( !M_linearSolver_u.converged() )
	  solveLinearSystemOnce_u( false );
      }

    if ( !M_linearSolver_u.converged() )
      std::cerr << "        WARNING: Velocity solver failed to converge."
		<< std::endl;

    return M_linearSolver_u.iterations();

} // solveLinearSystem



  template <typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::solveLinearSystemOnce_u( bool reusePC )
  {
    Chrono chrono;
    if ( reusePC )
      {
	chrono.start();
	M_linearSolver_u.solve( this->_u, M_f_u, SolverAztec::SAME_PRECONDITIONER );
	chrono.stop();

        if ( M_nUsePC_u == 1 )
	  {
	    M_tThisSolve_u = chrono.diff();
	    M_tLastSolve_u = M_tThisSolve_u;
	  }
	else
	  {
            M_tLastSolve_u = M_tThisSolve_u;
            M_tThisSolve_u = chrono.diff();
	  }
        M_tTotalSolve_u += M_tThisSolve_u;
        ++M_nUsePC_u;
      }
    else
      {
        chrono.start();
        M_linearSolver_u.solve( this->_u, M_f_u,
                              SolverAztec::SAME_NONZERO_PATTERN );
        chrono.stop();
        M_tThisSolve_u = chrono.diff();
        M_tLastSolve_u = 2 * M_tThisSolve_u;
        M_tTotalSolve_u = M_tThisSolve_u;
        M_nUsePC_u = 1;
      }
  } // solveLinearSystemOnce_u



  template <typename Mesh>
  UInt NavierStokesAleSolverCT<Mesh>::solveLinearSystem_p()
  {

    M_linearSolver_p.setMatrix( M_C_p );
    M_linearSolver_p.setRecursionLevel( 0 );

    if ( M_tTotalSolve_p < M_nUsePC_p * M_tThisSolve_p*M_tThisSolve_p/M_tLastSolve_p )
      solveLinearSystemOnce_p( false );
    else
      {
        solveLinearSystemOnce_p( true );
        if ( !M_linearSolver_p.converged() )
	  solveLinearSystemOnce_p( false );
      }

    if ( !M_linearSolver_p.converged() )
      std::cerr << "        WARNING: Pressure solver failed to converge."
		<< std::endl;

    return M_linearSolver_p.iterations();

} // solveLinearSystem



  template <typename Mesh>
  void NavierStokesAleSolverCT<Mesh>::solveLinearSystemOnce_p( bool reusePC )
  {
    Chrono chrono;
    if ( reusePC )
      {
	chrono.start();
	M_linearSolver_p.solve( this->_p, M_f_p, SolverAztec::SAME_PRECONDITIONER );
	chrono.stop();

        if ( M_nUsePC_p == 1 )
	  {
	    M_tThisSolve_p = chrono.diff();
	    M_tLastSolve_p = M_tThisSolve_p;
	  }
	else
	  {
            M_tLastSolve_p = M_tThisSolve_p;
            M_tThisSolve_p = chrono.diff();
	  }
        M_tTotalSolve_p += M_tThisSolve_p;
        ++M_nUsePC_p;
      }
    else
      {
        chrono.start();
        M_linearSolver_p.solve( this->_p, M_f_p,
                              SolverAztec::SAME_NONZERO_PATTERN );
        chrono.stop();
        M_tThisSolve_p = chrono.diff();
        M_tLastSolve_p = 2 * M_tThisSolve_p;
        M_tTotalSolve_p = M_tThisSolve_p;
        M_nUsePC_p = 1;
      }
  } // solveLinearSystemOnce_p

}
#endif
