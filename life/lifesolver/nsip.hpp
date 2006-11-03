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
  \file NavierStokesSolverIP.hpp
  \author M.A. Fernandez, Ch. Winkelmann
  \date 6/2003

  \brief This file contains a Navier-Stokes solver class which implements a
  stabilized implicit scheme with coupled solving.
*/
#ifndef _NAVIERSTOKESSOLVERIP_H_
#define _NAVIERSTOKESSOLVERIP_H_

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

#if AZTEC_SOLVER
#include <life/lifealg/SolverAztec.hpp>
#else

#include <life/lifearray/boostmatrix.hpp>
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
#include <life/lifesolver/nsipterms.hpp>
//#include <life/lifesolver/AFSolvers.hpp>

namespace LifeV
{
/*!
  \class NavierStokesSolverIP

  This class implements a NavierStokes solver via interior penalty
  stabilization. The resulting linear systems are solved by GMRES on the full
  matrix ( u and p coupled ).

*/
template<typename Mesh>
class NavierStokesSolverIP:
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
    NavierStokesSolverIP( const GetPot& dataFile,
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

    //! Initialize with steady Stokes solution (without convection)
    void initializeStokes( source_type const& source, Real t0 );

    //! linearize convective term around given (exact) velocity function
    void linearize( const Function& betaFct ) { M_betaFct = &betaFct; }

    //! removes mean of component comp of vector x
    void removeMean( Vector& x, UInt comp=1 );

    //! returns the residual vector
    Vector residual();

    //! calculates boundary force on a boundary specified by BCBase object
    void calculateBoundaryForce( EntityFlag flag,
                                 Real& fx , Real& fy, Real& fz );

private:

    /*! solve linear system
      @return number of iterations, zero if not applicable
      @param condEst condition estimate is returned here, -1 if not applicable
    */
    UInt solveLinearSystem( Real& condEst );

    //! solve linear system once (iterative solvers)
    void solveLinearSystemOnce( bool reusePC );

    //! apply boundary conditions
    void applyBoundaryConditions();

    //! Block pattern of M_u
    //MSRPatt M_pattMassUblock;

    //! Pattern for mass matrix
    //MixedPattern<nDimensions, nDimensions, MSRPatt> M_pattMassU;

    //! Block pattern of full matrix

#if AZTEC_SOLVER
    MSRPatt M_fullPattern;
    typedef MSRMatr<double> matrix_type;
#elif PETSC_SOLVER
    CSRPatt M_fullPattern;
    typedef BoostMatrix<boost::numeric::ublas::row_major> matrix_type;
#elif UMFPACK_SOLVER
    MSRPatt M_fullPattern;
    typedef BoostMatrix<boost::numeric::ublas::column_major> matrix_type;
#endif

    //! Mass matrix
    //MixedMatr<nDimensions, nDimensions, MSRPatt, double> M_matrMass;
    matrix_type M_matrMass;

    //! Stokes matrix: mu*stiff
    matrix_type M_matrStokes;

    //! Constant matrix: rho*bdfCoeff*mass + Stokes matrix
    matrix_type M_matrConst;

    //! Full matrix : constant matrix + convective term + stabilizations
    matrix_type M_matrFull;

    //! Full matrix without boundary conditions
    matrix_type M_matrNoBC;

    //! Elementary matrices and vectors
    ElemMat M_elmatStiff;
    ElemMat M_elmatBdfMass;
    ElemMat M_elmatMass;
    ElemMat M_elmatDiv;
    ElemMat M_elmatGrad;
    ElemVec M_elvec; // Elementary right hand side

    //! Right hand side for the velocity
    Vector M_rhsNoBC;

    //! Right hand side global
    Vector M_rhsFull;

    //! Global solution _u and _p
    Vector M_sol;

#if AZTEC_SOLVER
    typedef SolverAztec solver_type;
#elif PETSC_SOLVER
    typedef SolverPETSC solver_type;
#elif UMFPACK_SOLVER
    typedef SolverUMFPACK solver_type;
#endif
    solver_type M_linearSolver;

    Real M_time;

    bool M_steady;

    Real M_gammaBeta;
    Real M_gammaDiv;
    Real M_gammaPress;

    const Function* M_betaFct;

    bool M_divBetaUv;
    double M_diagonalize;

    Vector M_constantPressure;

    UInt M_nUsePC;
    Real M_tTotalSolve;
    Real M_tLastSolve;
    bool M_reusePC;

    UInt dim_u() const { return this->_dim_u; }
    UInt dim_p() const { return this->_dim_p; }
    BdfNS& bdf() { return this->_bdf; }
    const QuadRule& qr_u() const { return this->_Qr_u; }

}; // class NavierStokesSolverIP


template<typename Mesh> void NavierStokesSolverIP<Mesh>::
eval( Vector& fx0, Vector& /* gx0 */, Vector /* x0 */, int /* status */ )
{
    iterate( 0.0 );
    for ( UInt iDof = 0; iDof < nDimensions*dim_u() ; ++iDof )
    {
        fx0[ iDof ] = this->u()[ iDof ];
    }
}

//
// IMPLEMENTATION
//
template<typename Mesh> NavierStokesSolverIP<Mesh>::
NavierStokesSolverIP( const GetPot& dataFile,
                      const RefFE& refFE,
                      const QuadRule& quadRule,
                      const QuadRule& boundaryQuadRule,
                      BCHandler& bcHandler ):
    NavierStokesHandler<Mesh>( dataFile, refFE, refFE, quadRule,
                               boundaryQuadRule, quadRule, boundaryQuadRule,
                               bcHandler ),
    //M_pattMassUblock( uDof() ),
    //M_pattMassU( M_pattMassUblock, "diag" ),
#if AZTEC_SOLVER
    M_fullPattern( this->uDof(), this->mesh(), nDimensions+1 ),
#elif PETSC_SOLVER
    M_fullPattern( this->uDof(), nDimensions+1, this->mesh() ),
#elif UMFPACK_SOLVER
    M_fullPattern( this->uDof(), this->mesh(), nDimensions+1 ),
#endif
    //M_matrMassU( M_pattMassU ),
    M_matrMass( M_fullPattern ),
    M_matrStokes( M_fullPattern ),
    M_matrConst( M_fullPattern ),
    M_matrFull( M_fullPattern ),
    M_matrNoBC( M_fullPattern ),
    M_elmatStiff  ( this->fe_u().nbNode, nDimensions, nDimensions ),
    M_elmatBdfMass( this->fe_u().nbNode, nDimensions, nDimensions ),
    M_elmatMass   ( this->fe_u().nbNode, nDimensions, nDimensions ),
    M_elmatDiv    ( this->fe_u().nbNode, nDimensions+1, nDimensions ),
    M_elmatGrad   ( this->fe_u().nbNode, nDimensions, nDimensions+1 ),
    M_elvec( this->fe_u().nbNode, nDimensions ),
    M_rhsNoBC( ( nDimensions+1 )*dim_u() ),
    M_rhsFull( ( nDimensions+1 )*dim_u() ),
    M_sol( ( nDimensions+1 )*dim_u() ),
    M_betaFct( 0 ),
    M_constantPressure( ( nDimensions+1 )*dim_u() ),
    M_nUsePC( 1 ),
    M_tTotalSolve( 0 ),
    M_tLastSolve( 1 ),
    M_reusePC( false )
{
    M_steady = dataFile( "fluid/miscellaneous/steady", 1 );
    M_gammaBeta = dataFile( "fluid/ipstab/gammaBeta", 0. );
    M_gammaDiv = dataFile( "fluid/ipstab/gammaDiv", 0. );
    M_gammaPress = dataFile( "fluid/ipstab/gammaPress", 0. );
    M_divBetaUv = dataFile( "fluid/discretization/div_beta_u_v", 0);
    M_diagonalize = dataFile( "fluid/discretization/diagonalize", 1.);

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

#if AZTEC_SOLVER
    M_linearSolver.setOptionsFromGetPot( dataFile, "fluid/aztec" );
#elif PETSC_SOLVER
    M_linearSolver.setOptionsFromGetPot( dataFile, "fluid/petsc" );

    if ( this->bcHandler().hasOnlyEssential() && !M_diagonalize )
    {
        Real constPress = 1. / sqrt( dim_u() );
        for( UInt i=0; i<dim_u()*nDimensions; ++i )
        {
            M_constantPressure[ i ] = 0;
        }
        for( UInt i=dim_u()*nDimensions; i<dim_u()*(1+nDimensions); ++i )
        {
            M_constantPressure[ i ] = constPress;
        }
        std::vector<const Vector*> nullSpace(1);
        nullSpace[ 0 ] = &M_constantPressure;
        M_linearSolver.setNullSpace(nullSpace);
    }
#endif

    std::cout << std::endl;
    std::cout << "O-  Pressure unknowns: " << dim_p()     << std::endl;
    std::cout << "O-  Velocity unknowns: " << dim_u()     << std::endl
              << std::endl;
    std::cout << "O-  Computing constant matrices...        " << std::flush;

    Chrono chrono;
    chrono.start();

    // Matrices initialization
    //M_u.zeros();
    M_matrMass.zeros();
    M_matrStokes.zeros();
    M_matrConst.zeros();
    //chrono.stop();
    //std::cout << "(zeros:" << chrono.diff() << ") " << std::flush;
    //chrono.start();

    // Number of velocity components
    UInt nbCompU = this->u().nbcomp();

    Real bdfCoeff = bdf().bdf_u().coeff_der( 0 ) / this->dt();

    // Elementary computation and matrix assembling
    // Loop on elements
    for ( UInt iVol = 1; iVol <= this->mesh().numVolumes(); iVol++ )
    {
        this->fe_u().updateFirstDeriv( this->mesh().volumeList( iVol ) );

        M_elmatStiff.zero();
        M_elmatBdfMass.zero();
        M_elmatMass.zero();
        M_elmatDiv.zero();
        M_elmatGrad.zero();

        // stiffness strain computation
        stiff_strain( 2.0*this->viscosity(), M_elmatStiff, this->fe_u() );
        //stiff_div( 0.5*fe_u().diameter(), M_elmatStiff, fe_u() );

        // mass computation
        if ( !M_steady )
        {
            mass( this->density(), M_elmatMass, this->fe_u(), 0, 0, nDimensions );
            M_elmatBdfMass.mat() = M_elmatMass.mat();
            M_elmatBdfMass.mat() *= bdfCoeff;
        }

        for ( UInt iComp = 0; iComp<nbCompU; iComp++ )
        {
            for ( UInt jComp = 0; jComp<nbCompU; jComp++ )
            {
                // stiffness strain assembly
                assemb_mat( M_matrConst, M_elmatStiff, this->fe_u(),
                            this->uDof(), iComp, jComp );
                assemb_mat( M_matrStokes, M_elmatStiff, this->fe_u(),
                            this->uDof(), iComp, jComp );
            }

            // mass assembly
            if ( !M_steady )
            {
                assemb_mat( M_matrMass, M_elmatMass, this->fe_u(),
                            this->uDof(), iComp, iComp);
                assemb_mat( M_matrConst, M_elmatBdfMass, this->fe_u(),
                            this->uDof(), iComp, iComp);
            }

            // grad and div computation
            grad( iComp,  1.0, M_elmatGrad, this->fe_u(), this->fe_u(), iComp,
                  nDimensions);
            div ( iComp, -1.0, M_elmatDiv,  this->fe_u(), this->fe_u(), nbCompU,
                  iComp );

            // grad and div assembly
            assemb_mat( M_matrConst,  M_elmatGrad, this->fe_u(), this->uDof(),
                        iComp, nbCompU );
            assemb_mat( M_matrConst,  M_elmatDiv,  this->fe_u(), this->uDof(),
                        nbCompU, iComp );
            assemb_mat( M_matrStokes, M_elmatGrad, this->fe_u(), this->uDof(),
                        iComp, nbCompU );
            assemb_mat( M_matrStokes, M_elmatDiv,  this->fe_u(), this->uDof(),
                        nbCompU, iComp );

        }
    }

    //details::IPStabilization<Mesh, Dof>
    //    pressureStab(this->_mesh, this->uDof(), refFEu(), feBd_u(), qr_u(),
    //                 0, 0, M_gammaPress, this->viscosity() );
    //pressureStab.apply(M_matrConst, u());

    M_sol = ZeroVector( M_sol.size() );

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
    //M_matrStokes.spy( "stokes.m" );
    //M_matrConst.spy( "const.m" );
    //M_matrMass.spy("mass.m");


} // NavierStokesSolverIP() [Constructor]

template<typename Mesh>
void NavierStokesSolverIP<Mesh>::
timeAdvance( source_type const& source, Real const& time )
{
    M_time = time;

    std::cout << std::endl;
    std::cout << "O== Now we are at time "<< M_time << " s." << std::endl;

    // Number of velocity components
    UInt nbCompU = this->u().nbcomp();

    std::cout << "  o-  Updating mass term on right hand side... "
              << std::flush;

    Chrono chrono;
    chrono.start();

    // Right hand side for the velocity at time
    M_rhsNoBC = ZeroVector( M_rhsNoBC.size() );

    // loop on volumes: assembling source term
    for ( UInt iVol = 1; iVol<= this->mesh().numVolumes(); ++iVol )
    {
        M_elvec.zero();
        this->fe_u().updateJacQuadPt( this->mesh().volumeList( iVol ) );

        for ( UInt iComp = 0; iComp<nbCompU; ++iComp )
        {
            // compute local vector
            compute_vec( source, M_elvec, this->fe_u(), M_time, iComp );

            // assemble local vector into global one
            assemb_vec( M_rhsNoBC, M_elvec, this->fe_u(), this->uDof(), iComp );
        }
    }

    if ( !M_steady )
    {
        M_rhsNoBC += M_matrMass * bdf().bdf_u().time_der( this->dt() );
        //bdf().bdf_u().shift_right( M_sol );
    }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
} // timeAdvance()


template<typename Mesh>
void NavierStokesSolverIP<Mesh>::iterate( const Real& time )
{
    Chrono chrono;

    // velocity vector for linearization of convective term
    Vector betaVec( this->u().size() );

    if ( M_betaFct )
    {
        this->uInterpolate( *M_betaFct, betaVec, time );
    }
    else
    {
        if ( M_steady )
        {
            betaVec = this->u(); // last iteration
        }
        else
        {
            betaVec = bdf().bdf_u().extrap(); // bdf extrapolation
        }
    }

    M_time = time;

    std::cout << "  o-  Copying constant matrix...               "
              << std::flush;
    chrono.start();
    M_matrFull = M_matrConst;
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s."
              << std::endl;


    std::cout << "  o-  Updating convective volume terms...      "
              << std::flush;

    // Number of velocity components
    UInt nbCompU = this->u().nbcomp();

    chrono.start();

    // loop on volumes
    for ( UInt iVol = 1; iVol<= this->mesh().numVolumes(); ++iVol )
    {
        this->fe_u().updateFirstDeriv( this->mesh().volumeList( iVol ) ); //as updateFirstDer

        M_elmatStiff.zero();

        UInt eleID = this->fe_u().currentId();
        // Non linear term, Semi-implicit approach
        // M_elvec contains the velocity values in the nodes
        for ( UInt iNode = 0 ; iNode<( UInt )this->fe_u().nbNode ; iNode++ )
        {
            UInt  iloc = this->fe_u().patternFirst( iNode );
            for ( UInt iComp = 0; iComp<nbCompU; ++iComp )
            {
                UInt ig = this->uDof().localToGlobal( eleID, iloc+1 ) - 1 +
                    iComp*dim_u();
                M_elvec.vec()[ iloc+iComp*this->fe_u().nbNode ] =
                    this->density() * betaVec(ig);
            }
        }

        // Stabilising term: div u^n u v
        if ( M_divBetaUv )
        {
            mass_divw( 0.5*this->density(), M_elvec, M_elmatStiff, this->fe_u(),
                       0, 0, nbCompU );
        }

        // loop on components
        for ( UInt iComp = 0; iComp<nbCompU; ++iComp )
        {
            // compute local convective term and assembling
            grad( 0, M_elvec, M_elmatStiff, this->fe_u(), this->fe_u(),
                  iComp, iComp );
            grad( 1, M_elvec, M_elmatStiff, this->fe_u(), this->fe_u(),
                  iComp, iComp );
            grad( 2, M_elvec, M_elmatStiff, this->fe_u(), this->fe_u(),
                  iComp, iComp );
            assemb_mat( M_matrFull, M_elmatStiff, this->fe_u(), this->uDof(),
                        iComp, iComp);
        }
    }

    chrono.stop();

    std::cout << "done in " << chrono.diff() << " s."
              << std::endl;


    std::cout << "  o-  Updating convective ip terms...          "
              << std::flush;
    chrono.start();
    details::IPStabilization<Mesh, Dof>
        allStab( this->mesh(), this->uDof(), this->refFEu(), this->feBd_u(), qr_u(),
                 M_gammaBeta, M_gammaDiv, M_gammaPress, this->viscosity() );
    allStab.apply( M_matrFull, betaVec );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    M_matrNoBC = M_matrFull;

    // for BC treatment ( done at each time-step )
    std::cout << "  o-  Applying boundary conditions...          "
              << std::flush;
    chrono.start();
    applyBoundaryConditions();
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    // set initial guess
    if ( M_steady )
    {
        // use last iteration value as initial guess
        for ( UInt i = 0; i<nDimensions*dim_u(); ++i )
        {
            M_sol[ i ] = this->u()[ i ];
        }
    }
    else
    {
        // use bdf based extrapolation as initial guess
        M_sol = bdf().bdf_u().extrap();
    }
    if ( this->bcHandler().hasOnlyEssential() && !M_diagonalize )
    {
        removeMean( M_sol, 4 );
    }

    std::cout << "  o-  Solving system...                        "
              << std::flush;
    chrono.start();
    Real condEst;
    UInt iter = solveLinearSystem( condEst );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
    if ( iter > 0 )
    {
        std::cout << "        number of iterations                        = "
                  << iter << std::endl;

    }

    if ( condEst >= 0 )
    {
        std::cout << "        estimated condition number (preconditioned) = "
                  << condEst << std::endl;
    }

    if ( this->bcHandler().hasOnlyEssential() && !M_diagonalize )
    {
        removeMean( M_sol, 4 );
    }
    for ( UInt iDof = 0; iDof<nDimensions*dim_u(); ++iDof )
    {
        this->u()[ iDof ] = M_sol[ iDof ];
    }
    for ( UInt iDof = 0; iDof<dim_u(); ++iDof )
    {
        this->p()[ iDof ] = M_sol[ iDof+nDimensions*dim_u() ];
    }

    if ( !M_steady )
    {
        bdf().bdf_u().shift_right( M_sol );
    }

} // iterate()

template <typename Mesh>
void NavierStokesSolverIP<Mesh>::initialize( const Function& x0,
                                             Real t0, Real dt )
{
    ID nbComp = this->u().nbcomp(); // Number of components of the velocity
    bdf().bdf_u().initialize_unk( x0, this->mesh(), this->refFEu(), this->fe_u(), this->uDof(), t0,
                                 dt, nbComp+1 );

    // initialize M_sol with the first element in bdf_u.unk (=last value)
    M_sol = *( bdf().bdf_u().unk().begin() );
    if ( this->bcHandler().hasOnlyEssential() && !M_diagonalize )
    {
        removeMean( M_sol, 4 );
    }
    for ( UInt iDof = 0; iDof<nDimensions*dim_u(); ++iDof )
    {
        this->u()[ iDof ] = M_sol[ iDof ];
    }
    for ( UInt iDof = 0; iDof<dim_u(); ++iDof )
    {
        this->p()[ iDof ] = M_sol[ iDof+nDimensions*dim_u() ];
    }

    // write initial values (for debugging only)
    const std::vector<Vector>& unk = bdf().bdf_u().unk();
    for ( int iUnk=unk.size()-1; iUnk>=0; --iUnk )
    {
        for ( UInt iDof=0; iDof<nDimensions*dim_u(); ++iDof )
        {
            this->u()[ iDof ] = unk[ iUnk ][ iDof ];
        }
        for ( UInt iDof = 0; iDof<dim_u(); ++iDof )
        {
            this->p()[ iDof ] = unk[ iUnk ][ iDof+nDimensions*dim_u() ];
        }
        this->postProcess();
    }
} // initialize()

/*! Initialize when initial values for the velocity and the pressure are read
  from file
  @author Martin Prosi
*/
template <typename Mesh>
void
NavierStokesSolverIP<Mesh>::initialize( const std::string & vname )
{

    std::fstream resfile( vname.c_str(), std::ios::in | std::ios::binary );
    if ( resfile.fail() )
    {
        std::cerr << " Error in initialization: File not found or locked"
                  << std::endl;
        abort();
    }

    resfile.read( ( char* ) & M_sol( 0 ), M_sol.size() * sizeof( double ) );
    resfile.close();

    for ( UInt iDof = 0; iDof<nDimensions*dim_u(); ++iDof )
    {
        this->u()[ iDof ] = M_sol[ iDof ];
    }
    for ( UInt iDof = 0; iDof<dim_u(); ++iDof )
    {
        this->p()[ iDof ] = M_sol[ iDof+nDimensions*dim_u() ];
    }

    bdf().bdf_u().initialize_unk( M_sol );
    //bdf().bdf_p().initialize_unk( _p() );

    bdf().bdf_u().showMe();
    //bdf().bdf_p().showMe();

}

template <typename Mesh>
void
NavierStokesSolverIP<Mesh>::initializeStokes( source_type const& source,
                                              Real t0 )
{

    bool isSteady = M_steady;
    M_steady = true;
    timeAdvance( source, t0 );
    M_steady = isSteady;

    Chrono chrono;

    // velocity vector for linearization of convective term
    Vector betaVec(this->u().size());
    betaVec = ZeroVector( betaVec.size() );
    M_time = t0;

    std::cout << "  o-  Copying Stokes matrix...                 "
              << std::flush;
    chrono.start();

    M_matrFull = M_matrStokes;

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s."
              << std::endl;

    std::cout << "  o-  Adding ip terms...                       "
              << std::flush;
    chrono.start();
    details::IPStabilization<Mesh, Dof>
        initStab( this->mesh(), this->uDof(), this->refFEu(), this->feBd_u(), qr_u(),
                  0, 0, M_gammaPress, this->viscosity() );
    initStab.apply( M_matrFull, betaVec );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    std::cout << "  o-  Applying boundary conditions...          "
              << std::flush;
    chrono.start();
    applyBoundaryConditions();
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    M_sol = ZeroVector( M_sol.size() );

    std::cout << "  o-  Solving system...                        "
              << std::flush;
    chrono.start();
    Real condEst;
    UInt iter = solveLinearSystem( condEst );
    chrono.stop();

    // force rebuild next time
    M_reusePC = false;

    std::cout << "done in " << chrono.diff() << " s." << std::endl;
    if ( iter > 0 )
    {
        std::cout << "        number of iterations                        = "
                  << iter << std::endl;

    }
    if ( condEst >= 0 )
    {
        std::cout << "        estimated condition number (preconditioned) = "
                  << condEst << std::endl;
    }

    if ( this->bcHandler().hasOnlyEssential() && !M_diagonalize )
    {
        removeMean( M_sol, 4 );
    }
    for ( UInt iDof = 0; iDof<nDimensions*dim_u(); ++iDof )
    {
        this->u()[ iDof ] = M_sol[ iDof ];
    }
    for ( UInt iDof = 0; iDof<dim_u(); ++iDof )
    {
        this->p()[ iDof ] = M_sol[ iDof+nDimensions*dim_u() ];
    }

    bdf().bdf_u().initialize_unk( M_sol );

} // initializeStokes


template <typename Mesh>
void NavierStokesSolverIP<Mesh>::removeMean( Vector& x, UInt comp )
{
    Real sum1 = 0.;
    Real sum0 = 0.;
    for ( UInt iVol = 1; iVol <= this->mesh().numVolumes(); iVol++ )
    {
        this->fe_p().updateFirstDeriv( this->mesh().volumeList( iVol ) );
        sum1 += elem_integral( x, this->fe_p(), this->pDof(), comp );
        sum0 += this->fe_p().measure();
    }
    Real mean = sum1/sum0;
    for ( UInt iNode = (comp-1)*dim_u(); iNode<comp*dim_u(); ++iNode )
    {
        x[ iNode ] -= mean;
    }

} // removeMean()

template <typename Mesh>
void NavierStokesSolverIP<Mesh>::applyBoundaryConditions()
{
    //     M_rhsFull = 0.0;
    //     for ( UInt i = 0; i<nDimensions*dim_u(); ++i )
    //     {
    //         M_rhsFull[ i ] = M_rhsNoBC[ i ];
    //     }
    M_rhsFull = M_rhsNoBC;

    // BC manage for the velocity
    if ( !this->bcHandler().bdUpdateDone() )
        this->bcHandler().bdUpdate( this->mesh(), this->feBd_u(), this->uDof() );
    bcManage( M_matrFull, M_rhsFull, this->mesh(), this->uDof(), this->bcHandler(), this->feBd_u(), 1.0,
               M_time );

    if ( this->bcHandler().hasOnlyEssential() && M_diagonalize )
    {
         M_matrFull.diagonalize( nDimensions*dim_u(), M_diagonalize,
                                 M_rhsFull, 0);
        //                         pexact( M_time,
        //                                 this->_mesh.point( 1 ).x(),
        //                                 this->_mesh.point( 1 ).y(),
        //                                 this->_mesh.point( 1 ).z(), 1 ) );
        // M_matrFull.diagonalize_row( nDimensions*dim_u(), 1.0 );
        // M_rhsFull[ nDimensions*dim_u() ] = pexact( this->_mesh.point( 1 ).x(),
        //                                         this->_mesh.point( 1 ).y(),
        //                                         this->_mesh.point( 1 ).z() );
    }
} // applyBoundaryCondition

template <typename Mesh>
UInt NavierStokesSolverIP<Mesh>::solveLinearSystem(Real& condEst)
{
#if 1-UMFPACK_SOLVER
    M_linearSolver.setMatrix( M_matrFull );

    bool reuse = M_reusePC;
    solveLinearSystemOnce( M_reusePC );
    if ( reuse && !M_linearSolver.converged() ) // retry if not converged
        solveLinearSystemOnce( false );
#else
    M_linearSolver.setMatrix( M_matrFull, true );

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
void NavierStokesSolverIP<Mesh>::solveLinearSystemOnce( bool reusePC )
{
    Chrono chrono;
    if ( reusePC ) {
        chrono.start();
#if AZTEC_SOLVER
        M_linearSolver.solve( M_sol, M_rhsFull,
                              SolverAztec::SAME_PRECONDITIONER );
#elif PETSC_SOLVER
        M_linearSolver.solve( M_sol, M_rhsFull, SAME_PRECONDITIONER );
#endif
        chrono.stop();
        Real tThisSolve = chrono.diff();
        M_tTotalSolve += tThisSolve;
        ++M_nUsePC;

        Real tNextSolve;
        if ( M_nUsePC == 2 ) {
            tNextSolve = tThisSolve;
        } else {
            tNextSolve = tThisSolve * tThisSolve / M_tLastSolve;
        }
        M_reusePC = ( M_tTotalSolve > M_nUsePC * tNextSolve );
        M_tLastSolve = tThisSolve;

    } else {
        chrono.start();
#if AZTEC_SOLVER
        M_linearSolver.solve( M_sol, M_rhsFull,
                              SolverAztec::SAME_NONZERO_PATTERN );
#elif PETSC_SOLVER
        M_linearSolver.solve( M_sol, M_rhsFull, SAME_NONZERO_PATTERN );
#endif
        chrono.stop();
        M_tLastSolve = chrono.diff();
        M_reusePC = true;
        M_tTotalSolve = M_tLastSolve;
        M_nUsePC = 1;
    }
} // solveLinearSystemOnce

template <typename Mesh>
Vector NavierStokesSolverIP<Mesh>::residual()
{
    Vector residual( ( nDimensions+1 )*dim_u() );
    residual = M_rhsNoBC;
    residual -= M_matrNoBC * M_sol;
    return residual;
}

template <typename Mesh>
void NavierStokesSolverIP<Mesh>::calculateBoundaryForce( EntityFlag flag,
                                                         Real& fx,
                                                         Real& fy,
                                                         Real& fz )
{
    Vector residual( ( nDimensions+1 )*dim_u() );
    residual = M_matrNoBC * M_sol - M_rhsNoBC;

    //BCFunctionBase f( fct );
    //BCHandler bch;
    //bch.addBC( "force boundary x", flag, Essential, Full, f, nDimensions );
    //bch.bdUpdate( this->_mesh, feBd_u(), this->uDof() );
    //const BCBase& BCb = bch[0]

    const BCBase& BCb = this->bcHandler().GetBCWithFlag( flag );

    fx = 0;
    fy = 0;
    fz = 0;

    // Number of total scalar Dof
    const UInt nDof = this->uDof().numTotalDof();

    // Loop on BC identifiers
    for ( ID i = 1; i <= BCb.list_size(); ++i )
    {
        ID idDof = BCb( i )->id();

        fx += residual( idDof );
        fy += residual( idDof + nDof );
        fz += residual( idDof + 2*nDof );
    }

//     const UInt nDof = M_dof.numTotalDof();

//     // local trace of the residual
//     ElemVec f( M_feBd.nbNode, nDimensions );

//     // loop on boundary faces
//     for ( UInt iFace = 1; iFace<=M_mesh.numBFaces(); ++iFace )
//     {
//         // update current finite elements
//         M_feBd.updateMeas( M_mesh.face( iFace ) );
//         const Real area = M_feBd.measure();
//         const UInt iElAd1 = M_mesh.face( iFace ).ad_first();
//         M_fe1.updateFirstDeriv( M_mesh.volumeList( iElAd1 ) );

//         // local id of the face in its adjacent element
//         UInt iFaEl = M_mesh.face( iFace ).pos_first();
//         for ( int iNode = 0; iNode < M_feBd.nbNode; ++iNode )
//         {
//             UInt iloc = M_fToP( iFaEl, iNode+1 );
//             for ( int iCoor = 0; iCoor < M_fe1.nbCoor; ++iCoor )
//             {
//                 UInt ig = this->uDof().localToGlobal( iElAd1, iloc+1 )-1+iCoor*nDof;
//                 f.vec()[ iCoor*M_feBd.nbNode + iNode ] = state( ig );
//             }
//         }

//     } // loop on boundary faces

} // calculateBoundaryForce

} // namespace LifeV

#endif //_NAVIERSTOKESSOLVERIP_H_
