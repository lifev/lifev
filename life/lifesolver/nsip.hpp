/*!
  \file NavierStokesSolverIP.h
  \author M.A. Fernandez
  \date 6/2003
  \version 1.0

  \brief This file contains a Navier-Stokes solver class which implements a
  stabilized implicit scheme with coupled solving.
*/
#ifndef _NAVIERSTOKESSOLVERIP_H_
#define _NAVIERSTOKESSOLVERIP_H_

#define USE_AZTEC_SOLVER 1

#include <life/lifesolver/NavierStokesHandler.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>

#if USE_AZTEC_SOLVER
#include <life/lifealg/SolverAztec.hpp>
#else
#include <lifeconfig.h>
#if defined( HAVE_PETSC_H )
#include <life/lifealg/SolverPETSC.hpp>
#endif /* HAVE_PETSC_H */
#endif /* USE_AZTEC_SOLVER */

#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifefem/sobolevNorms.hpp>
#include <life/lifefem/geoMap.hpp>
#include <ipStabilization.hpp>

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

    //! linearize convective term around given (exact) velocity function
    void linearize( const Function& betaFct ) { M_betaFct = &betaFct; }

    //! removes mean of component comp of vector x
    void removeMean( Vector& x, UInt comp=1 );

private:

    //! Block pattern of M_u
    //MSRPatt M_pattMassUblock;

    //! Pattern for mass matrix
    //MixedPattern<nDimensions, nDimensions, MSRPatt> M_pattMassU;

    //! Block pattern of full matrix
    MSRPatt M_fullPattern;

    //! Matrix M_u: Vmass
    //MixedMatr<nDimensions, nDimensions, MSRPatt, double> M_matrMass;
    MSRMatr<double> M_matrMass;

    //! Matrix CStokes: rho*bdfCoeff*Vmass + mu*Vstiff + linear stabilizations
    MSRMatr<double> M_matrStokes;

    //! Matrix C: CStokes + Convective_term + nonlinear stabilizations
    MSRMatr<double> M_matrFull;

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
    Vector M_rhsFull;

    //! Global solution _u and _p
    Vector M_sol;

#if USE_AZTEC_SOLVER
    SolverAztec M_linearSolver;
#else
    SolverPETSC M_linearSolver;
#endif

    Real M_time;

    bool M_steady;

    Real M_gammaBeta;
    Real M_gammaDiv;
    Real M_gammaPress;

    const Function* M_betaFct;

    bool M_divBetaUv;
    double M_diagonalize;

    Vector M_constantPressure;

}; // class NavierStokesSolverIP


template<typename Mesh> void NavierStokesSolverIP<Mesh>::
eval( Vector& fx0, Vector& /* gx0 */, Vector /* x0 */, int /* status */ )
{
    iterate( 0.0 );
    for ( UInt iDof = 0; iDof < nDimensions*_dim_u ; ++iDof )
    {
        fx0[ iDof ] = _u[ iDof ];
    }
}

//
//                                         IMPLEMENTATION
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
    //M_pattMassUblock( _dof_u ),
    //M_pattMassU( M_pattMassUblock, "diag" ),
    M_fullPattern( _dof_u, _mesh, nDimensions+1 ),
    //M_matrMassU( M_pattMassU ),
    M_matrMass( M_fullPattern ),
    M_matrStokes( M_fullPattern ),
    M_matrFull( M_fullPattern ),
    M_elmatC( _fe_u.nbNode, nDimensions, nDimensions ),
    M_elmatMass( _fe_u.nbNode, nDimensions, nDimensions ),
    M_elmatD( _fe_u.nbNode, nDimensions+1, nDimensions ),
    M_elmatDtr( _fe_u.nbNode, nDimensions, nDimensions+1 ),
    M_elmatP( _fe_u.nbNode, nDimensions+1, nDimensions+1 ),
    M_elvec( _fe_u.nbNode, nDimensions ),
    M_rhsU( ( nDimensions+1 )*_dim_u ),
    M_rhsFull( ( nDimensions+1 )*_dim_u ),
    M_sol( ( nDimensions+1 )*_dim_u ),
    M_betaFct( 0 ),
    M_constantPressure( ( nDimensions+1 )*_dim_u )
{
    M_steady = dataFile( "fluid/miscellaneous/steady", 1 );
    M_gammaBeta = dataFile( "fluid/ipstab/gammaBeta", 0. );
    M_gammaDiv = dataFile( "fluid/ipstab/gammaDiv", 0. );
    M_gammaPress = dataFile( "fluid/ipstab/gammaPress", 0. );
    M_divBetaUv = dataFile( "fluid/discretization/div_beta_u_v", 0);
    M_diagonalize = dataFile( "fluid/discretization/diagonalize", 1.);

    // check mesh for elements with all nodes on the boundary
    UInt nLocalFaces = _mesh.numLocalFaces();
    UInt nFixedTets = 0;
    for ( UInt iVol = 1; iVol <= _mesh.numVolumes(); iVol++ )
    {
        UInt nBoundaryFaces = 0;
        for ( UInt iLocalFace=1; iLocalFace<=nLocalFaces; ++iLocalFace )
        {
            if ( _mesh.isBoundaryFace( _mesh.localFaceId( iVol, iLocalFace )))
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
        std::cerr << "WARNING: " << nFixedTets << " of " << _mesh.numVolumes()
                  << " tetrahedrons have all nodes on the boundary."
                  << std::endl;
    }

#if USE_AZTEC_SOLVER
    M_linearSolver.setOptionsFromGetPot( dataFile, "fluid/aztec" );
#else
    M_linearSolver.setOptionsFromGetPot( dataFile, "fluid/petsc" );

    if ( _BCh_u.hasOnlyEssential() && !M_diagonalize )
    {
        Real constPress = 1. / sqrt( _dim_u );
        for( UInt i=0; i<_dim_u*nDimensions; ++i )
        {
            M_constantPressure[ i ] = 0;
        }
        for( UInt i=_dim_u*nDimensions; i<_dim_u*(1+nDimensions); ++i )
        {
            M_constantPressure[ i ] = constPress;
        }
        std::vector<const Vector*> nullSpace(1);
        nullSpace[ 0 ] = &M_constantPressure;
        M_linearSolver.setNullSpace(nullSpace);
    }
#endif

    std::cout << std::endl;
    std::cout << "O-  Pressure unknowns: " << _dim_p     << std::endl;
    std::cout << "O-  Velocity unknowns: " << _dim_u     << std::endl
              << std::endl;
    std::cout << "O-  Computing mass and Stokes matrices... " << std::flush;

    Chrono chrono;
    chrono.start();

    // Matrices initialization
    //M_u.zeros();
    M_matrMass.zeros();
    M_matrStokes.zeros();

    // Number of velocity components
    UInt nbCompU = _u.nbcomp();

    Real bdfCoeff = _bdf.bdf_u().coeff_der( 0 ) / _dt;

    // Elementary computation and matrix assembling
    // Loop on elements
    for ( UInt iVol = 1; iVol <= _mesh.numVolumes(); iVol++ )
    {
        _fe_u.updateFirstDeriv( _mesh.volumeList( iVol ) );

        M_elmatC.zero();
        M_elmatMass.zero();
        M_elmatD.zero();
        M_elmatDtr.zero();

        // stiffness strain
        stiff_strain( 2.0*_mu, M_elmatC, _fe_u );
        //stiff_div( 0.5*_fe_u.diameter(), M_elmatC, _fe_u );

        // mass
        if ( !M_steady )
        {
            mass( _rho*bdfCoeff, M_elmatMass, _fe_u, 0, 0, nDimensions );
            M_elmatC.mat() += M_elmatMass.mat();
            M_elmatMass.mat() *= ( 1./bdfCoeff );
        }

        for ( UInt iComp = 0;iComp<nbCompU;iComp++ )
        {
            for ( UInt jComp = 0;jComp<nbCompU;jComp++ )
            {
                // stiffness
                assemb_mat( M_matrStokes, M_elmatC, _fe_u, _dof_u, iComp,
                            jComp );
                if ( !M_steady )
                {
                    assemb_mat( M_matrMass, M_elmatMass, _fe_u, _dof_u, iComp,
                                jComp);
                }
            }
            // mass
            //assemb_mat( M_matrMassU, M_elmatMass, _fe_u,_dof_u,iComp,iComp );

            // div
            grad( iComp, 1.0, M_elmatDtr, _fe_u, _fe_u, iComp, nDimensions );
            div( iComp, -1.0, M_elmatD  , _fe_u, _fe_u, nbCompU, iComp );

            // assembling
            assemb_mat( M_matrStokes, M_elmatDtr, _fe_u, _dof_u, iComp,
                        nbCompU );
            assemb_mat( M_matrStokes, M_elmatD, _fe_u, _dof_u, nbCompU,
                        iComp );

        }
    }

    //IPStabilization<Mesh, Dof>
    //    pressureStab(_mesh, _dof_u, _refFE_u, _feBd_u, _Qr_u,
    //                 0, 0, M_gammaPress, this->viscosity() );
    //pressureStab.apply(M_matrStokes, this->_u);

    M_sol = ZeroVector( M_sol.size() );

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
    // M_matrStokes.spy( "CS.m" );


} // NavierStokesSolverIP() [Constructor]

template<typename Mesh>
void NavierStokesSolverIP<Mesh>::
timeAdvance( source_type const& source, Real const& time )
{
    M_time = time;

    std::cout << std::endl;
    std::cout << "O== Now we are at time "<< M_time << " s." << std::endl;

    // Number of velocity components
    UInt nbCompU = _u.nbcomp();

    std::cout << "  o-  Updating mass term on right hand side... "
              << std::flush;

    Chrono chrono;
    chrono.start();

    // Right hand side for the velocity at time
    M_rhsU = ZeroVector( M_rhsU.size() );

    // loop on volumes: assembling source term
    for ( UInt iVol = 1; iVol<= _mesh.numVolumes(); ++iVol )
    {
        M_elvec.zero();
        _fe_u.updateJacQuadPt( _mesh.volumeList( iVol ) );

        for ( UInt iComp = 0; iComp<nbCompU; ++iComp )
        {
            // compute local vector
            compute_vec( source, M_elvec, _fe_u, M_time, iComp );

            // assemble local vector into global one
            assemb_vec( M_rhsU, M_elvec, _fe_u, _dof_u, iComp );
        }
    }

    if ( !M_steady )
    {
        M_rhsU += M_matrMass * _bdf.bdf_u().time_der( _dt );
        //_bdf.bdf_u().shift_right( M_sol );
    }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
} // timeAdvance()


template<typename Mesh>
void NavierStokesSolverIP<Mesh>::iterate( const Real& time )
{
    Chrono chrono;

    // velocity vector for linearization of convective term
    Vector betaVec(_u.size());

    if ( M_betaFct )
    {
        this->uInterpolate( *M_betaFct, betaVec, time );
    }
    else
    {
        if ( M_steady )
        {
            betaVec = _u; // last iteration
        }
        else
        {
            betaVec = _bdf.bdf_u().extrap(); // bdf extrapolation
        }
    }

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
    UInt nbCompU = _u.nbcomp();

    chrono.start();

    // loop on volumes
    for ( UInt iVol = 1; iVol<= _mesh.numVolumes(); ++iVol )
    {
        _fe_u.updateFirstDeriv( _mesh.volumeList( iVol ) ); //as updateFirstDer

        M_elmatC.zero();

        UInt eleID = _fe_u.currentId();
        // Non linear term, Semi-implicit approach
        // M_elvec contains the velocity values in the nodes
        for ( UInt iNode = 0 ; iNode<( UInt )_fe_u.nbNode ; iNode++ )
        {
            UInt  iloc = _fe_u.patternFirst( iNode );
            for ( UInt iComp = 0; iComp<nbCompU; ++iComp )
            {
                UInt ig = _dof_u.localToGlobal( eleID, iloc+1 )-1+iComp*_dim_u;
                M_elvec.vec()[ iloc+iComp*_fe_u.nbNode ] = _rho * betaVec(ig);
            }
        }

        // Stabilising term: div u^n u v
        if ( M_divBetaUv )
        {
            mass_divw( 0.5*_rho, M_elvec, M_elmatC, _fe_u, 0, 0, nbCompU );
        }

        // loop on components
        for ( UInt iComp = 0; iComp<nbCompU; ++iComp )
        {
            // compute local convective term and assembling
            grad( 0, M_elvec, M_elmatC, _fe_u, _fe_u, iComp, iComp );
            grad( 1, M_elvec, M_elmatC, _fe_u, _fe_u, iComp, iComp );
            grad( 2, M_elvec, M_elmatC, _fe_u, _fe_u, iComp, iComp );
            assemb_mat( M_matrFull, M_elmatC, _fe_u, _dof_u, iComp, iComp );
        }
    }

    chrono.stop();

    std::cout << "done in " << chrono.diff() << " s."
              << std::endl;

    std::cout << "  o-  Updating convective ip terms...          "
              << std::flush;

    chrono.start();

    IPStabilization<Mesh, Dof>
        allStab( _mesh, _dof_u, _refFE_u, _feBd_u, _Qr_u,
                 M_gammaBeta, M_gammaDiv, M_gammaPress, this->viscosity() );
    allStab.apply( M_matrFull, betaVec );

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;



    // for BC treatment ( done at each time-step )
    std::cout << "  o-  Applying boundary conditions...          "
              << std::flush;
    chrono.start();



    //     M_rhsFull = 0.0;
    //     for ( UInt i = 0; i<nDimensions*_dim_u; ++i )
    //     {
    //         M_rhsFull[ i ] = M_rhsU[ i ];
    //     }
    M_rhsFull = M_rhsU;

    // BC manage for the velocity
    if ( !_BCh_u.bdUpdateDone() )
        _BCh_u.bdUpdate( _mesh, _feBd_u, _dof_u );
    bcManage( M_matrFull, M_rhsFull, _mesh, _dof_u, _BCh_u, _feBd_u, 1.0,
               M_time );

    if ( _BCh_u.hasOnlyEssential() && M_diagonalize )
    {
         M_matrFull.diagonalize( nDimensions*_dim_u, M_diagonalize,
                                 M_rhsFull, 0);
        //                         pexact( M_time,
        //                                 _mesh.point( 1 ).x(),
        //                                 _mesh.point( 1 ).y(),
        //                                 _mesh.point( 1 ).z(), 1 ) );
        // M_matrFull.diagonalize_row( nDimensions*_dim_u, 1.0 );
        // M_rhsFull[ nDimensions*_dim_u ] = pexact( _mesh.point( 1 ).x(),
        //                                         _mesh.point( 1 ).y(),
        //                                         _mesh.point( 1 ).z() );
    }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    M_linearSolver.setMatrix( M_matrFull );



// ---------------
// C * V = F
// ---------------

    // set initial guess
    if ( M_steady )
    {
        // use last iteration value as initial guess
        for ( UInt i = 0; i<nDimensions*_dim_u; ++i )
        {
            M_sol[ i ] = _u[ i ];
        }
    }
    else
    {
        // use bdf based extrapolation as initial guess
        M_sol = _bdf.bdf_u().extrap();
    }
    removeMean( M_sol, 4 );

    std::cout << "  o-  Solving system...                        "
              << std::flush;
    chrono.start();
#if USE_AZTEC_SOLVER
    M_linearSolver.solve( M_sol, M_rhsFull,
                          SolverAztec::SAME_PRECONDITIONER );
#else
    M_linearSolver.solve( M_sol, M_rhsFull, SAME_PRECONDITIONER );
#endif
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    if ( !M_linearSolver.converged() )
    {
        std::cerr << "        WARNING: Solver failed to converge."
                  << std::endl;
    }
#if USE_AZTEC_SOLVER
#else
    std::cout << "        estimated condition number (preconditioned) = "
              << M_linearSolver.condEst() << std::endl;
#endif
    std::cout << "        number of iterations                        = "
              << M_linearSolver.iterations() << std::endl;

    removeMean( M_sol, 4 );
    for ( UInt iDof = 0; iDof<nDimensions*_dim_u; ++iDof )
    {
        _u[ iDof ] = M_sol[ iDof ];
    }
    for ( UInt iDof = 0; iDof<_dim_u; ++iDof )
    {
        _p[ iDof ] = M_sol[ iDof+nDimensions*_dim_u ];
    }

    if ( !M_steady )
    {
        _bdf.bdf_u().shift_right( M_sol );
    }

} // iterate()

template <typename Mesh>
void NavierStokesSolverIP<Mesh>::initialize( const Function& x0,
                                             Real t0, Real dt )
{
    ID nbComp = _u.nbcomp(); // Number of components of the velocity
    _bdf.bdf_u().initialize_unk( x0, this->_mesh, _refFE_u, _fe_u, _dof_u, t0,
                                 dt, nbComp+1 );

    // initialize M_sol with the first element in bdf_u.unk (=last value)
    M_sol = *( _bdf.bdf_u().unk().begin() );
    removeMean( M_sol, 4 );
    for ( UInt iDof = 0; iDof<nDimensions*_dim_u; ++iDof )
    {
        _u[ iDof ] = M_sol[ iDof ];
    }
    for ( UInt iDof = 0; iDof<_dim_u; ++iDof )
    {
        _p[ iDof ] = M_sol[ iDof+nDimensions*_dim_u ];
    }

    // write initial values (for debugging only)
    const std::vector<Vector>& unk = _bdf.bdf_u().unk();
    for ( int iUnk=unk.size()-1; iUnk>=0; --iUnk )
    {
        for ( UInt iDof=0; iDof<nDimensions*_dim_u; ++iDof )
        {
            _u[ iDof ] = unk[ iUnk ][ iDof ];
        }
        for ( UInt iDof = 0; iDof<_dim_u; ++iDof )
        {
            _p[ iDof ] = unk[ iUnk ][ iDof+nDimensions*_dim_u ];
        }
        this->postProcess();
    }
} // initialize()

template <typename Mesh>
void NavierStokesSolverIP<Mesh>::removeMean( Vector& x, UInt comp )
{
    Real sum1 = 0.;
    Real sum0 = 0.;
    for ( UInt iVol = 1; iVol <= _mesh.numVolumes(); iVol++ )
    {
        _fe_p.updateFirstDeriv( _mesh.volumeList( iVol ) );
        sum1 += elem_integral( x, _fe_p, _dof_p, comp );
        sum0 += _fe_p.measure();
    }
    Real mean = sum1/sum0;
    for ( UInt iNode = (comp-1)*_dim_u; iNode<comp*_dim_u; ++iNode )
    {
        x[ iNode ] -= mean;
    }

} // removeMean()

} // namespace LifeV

#endif //_NAVIERSTOKESSOLVERIP_H_
