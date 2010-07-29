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
  \author G. Fourestey
  \date 2/2007

  \brief This file contains a Navier-Stokes solver class
*/
#ifndef _ADRSOLVER_H_
#define _ADRSOLVER_H_

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>

#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>
//
#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifefem/sobolevNorms.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifesolver/nsipterms.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
//
#include <boost/shared_ptr.hpp>

#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>


namespace LifeV
{
/*!
  \class ADRSolver

  This class implements a NavierStokes solver.
  The resulting linear systems are solved by GMRES on the full
  matrix ( u and p coupled ).

*/




template< typename Mesh,
          typename SolverType = LifeV::Epetra::SolverTrilinos >
class ADRSolver
//     :
//     public NavierStokesHandler<Mesh>, EpetraHandler
{

public:

    typedef DataNavierStokes<Mesh> data_type;
    typedef Real ( *Function ) ( const Real&, const Real&, const Real&,
                                 const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> source_type;

    typedef Mesh mesh_type;

    typedef BCHandler                             bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type> bchandler_type;

//    typedef SolverType                    solver_type;
    typedef typename SolverType::matrix_type      matrix_type;
    typedef typename SolverType::vector_type      vector_type;



    //! Constructor
    /*!
      \param dataType
      \param localMap
      \param velocity FE space
      \param pressure FE space
      \param bcHandler boundary conditions for the velocity
    */
    ADRSolver( const data_type&     dataType,
               const EpetraMap&     localMap,
               FESpace<Mesh>&       uFESpace,
               FESpace<Mesh>&       pFESpace,
               BCHandler&           bcHandler,
               Epetra_Comm&         comm );

    //! virtual destructor

    virtual ~ADRSolver();

    //! Update the right hand side for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    virtual void timeAdvance( source_type const& source, Real const& time );

    //! Update convective term, bc treatment and solve the linearized ns system
    virtual void iterate( const Real& time );

    void eval( Vector& fx0, Vector& gx0, Vector x0, int status );

    void initialize( const Function& x0, Real t0=0., Real dt=0. );

    /*! Initialize when initial values for the velocity and the pressure are
      read from file
      @author Martin Prosi
    */
    void initialize( const std::string & vname );
    void initializeNull();

    //! Initialize with steady Stokes solution (without convection)
    void initializeStokes( source_type const& source, Real t0 );

//     //! linearize convective term around given (exact) velocity function
//     void linearize( const Function& betaFct ) { std::cout << "linearize " << std::endl;M_betaFct = &betaFct; }

    void setUp        ( const GetPot& dataFile );

    void buildSystem();
    void updateSystem();

    //! Returns the velocity vector
    PhysVectUnknown<Vector>& u() {return M_u;}

    //! Returns the pressure
    ScalUnknown<Vector>& p()     {return M_p;}

    //! Bounday Conditions
    const bool setFluidBC() const {return M_setBC;}
    //! set the fluid BCs
    void setFluidBC(BCHandler &BCh_u)
        {
            M_BCh_fluid = &BCh_u; M_setBC = true;
        }

    BCHandler& bcHandler()
        {
            return *M_BCh_fluid;
        }

    //! set the source term functor
    void setSourceTerm( source_type __s )
        {
            M_source = __s;
        }

    //! get the source term functor
    source_type sourceTerm() const
        {
            return M_source;
        }

    //! Postprocessing
    virtual void postProcess();

    //
    UInt dim_u() const           { return M_uFESpace.dim(); }
    UInt dim_p() const           { return M_pFESpace.dim(); }


protected:


    void solveSystem            ( );
    void applyBoundaryConditions( matrix_type& _csr  );

    void echo(std::string message);

    //! removes mean of component comp of vector x
    void removeMean( vector_type& x, UInt comp = 1 );



    //private members

    //! data for NS solvers
    const data_type&               M_dataType;

    // FE spaces
    FESpace<Mesh>&                 M_uFESpace;
    FESpace<Mesh>&                 M_pFESpace;

    //! bdf stuff
    BdfTNS<vector_type>            M_bdf;

    //! MPI communicator
    Epetra_Comm*                   M_comm;
    int                            M_me;

    //! fluid BC
    BCHandler*                     M_BCh_fluid;
    bool                           M_setBC;

    EpetraMap                      M_localMap;

    //! mass matrix

    matrix_type*                   M_matrMass;

    //! Stokes matrix: mu*stiff
    matrix_type*                   M_matrStokes;

    //! matrix to be solved
    matrix_type*                   M_matrFull;

    //! source term for NS
    source_type                    M_source;


    //! Right hand side for the velocity
    vector_type                    M_rhsNoBC;

    //! Right hand side global
    vector_type                    M_rhsFull;

    //! Global solution _u and _p
    vector_type                    M_sol;


    SolverType                     M_linearSolver;


    bool                           M_steady;

    //! Stabilization
    bool                           M_stab;

    details::
    IPStabilization<Mesh, Dof>     M_ipStab;
    Real                           M_gammaBeta;
    Real                           M_gammaDiv;
    Real                           M_gammaPress;

    const Function*                M_betaFct;

    bool                           M_divBetaUv;

    //
    double                         M_diagonalize;

    UInt                           M_count;

    //! boolean that indicates if output is sent to cout

    bool                           M_verbose;

    //! boolean that indicates if the matrix is updated for the current iteration

    bool                           M_updated;

    //! boolean that indicates if le precond has to be recomputed

    bool                           M_reusePrec;

    //! interger storing the max number of solver iteration with prec recomputing

    int                            M_maxIterSolver;

    //!

    bool                           M_recomputeMatrix;

private:

    //! Elementary matrices and vectors
    ElemMat                        M_elmatStiff;      // velocity Stokes
    ElemMat                        M_elmatBdfMass;    // velocity Stokes + coeff*mass
    ElemMat                        M_elmatMass;       // velocity mass
    ElemMat                        M_elmatP;          // (p,q) bloc for preconditioners
    ElemMat                        M_elmatDiv;
    ElemMat                        M_elmatGrad;
    ElemVec                        M_elvec;           // Elementary right hand side

    //! The velocity
    PhysVectUnknown<Vector>        M_u;

    //! The pressure
    ScalUnknown<Vector>            M_p;


}; // class ADRSolver


template<typename Mesh, typename SolverType>
void
ADRSolver<Mesh, SolverType>::eval( Vector& fx0, Vector& /* gx0 */, Vector /* x0 */, int /* status */ )
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


template<typename Mesh, typename SolverType>
ADRSolver<Mesh, SolverType>::
ADRSolver( const data_type&     dataType,
                    const EpetraMap&     localMap,
                    FESpace<Mesh>&       uFESpace,
                    FESpace<Mesh>&       pFESpace,
                    BCHandler&           BCh_u,
                    Epetra_Comm&         comm ):
    M_dataType               ( dataType ),
    M_uFESpace               ( uFESpace ),
    M_pFESpace               ( pFESpace ),
    M_bdf                    ( M_dataType.dataTime()->getBDF_order() ),
    M_BCh_fluid              ( &BCh_u ),
    M_comm                   ( &comm ),
    M_me                     ( comm.MyPID() ),
    M_linearSolver           ( ),
    M_localMap               ( localMap ),
    M_matrMass               ( 0 ),
    M_matrStokes             ( 0 ),
    M_matrFull               ( 0 ),
    M_elmatStiff             ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatBdfMass           ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatMass              ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatP                 ( M_pFESpace.fe().nbNode, 1, 1 ),
    M_elmatDiv               ( M_pFESpace.fe().nbNode, 1, 0, M_uFESpace.fe().nbNode, 0, nDimensions ),
    M_elmatGrad              ( M_uFESpace.fe().nbNode, nDimensions, 0, M_pFESpace.fe().nbNode, 0, 1 ),
    M_elvec                  ( M_uFESpace.fe().nbNode, nDimensions ),
    M_rhsNoBC                ( localMap ),
    M_rhsFull                ( localMap ),
    M_sol                    ( localMap ),
    M_stab                   ( false ),
    M_ipStab                 ( M_dataType.dataMesh()->mesh(),
                               M_uFESpace.dof(), M_uFESpace.refFE(),
                               M_uFESpace.feBd(), M_uFESpace.qr(),
                               0., 0., 0.,
                               M_dataType.viscosity() ),
    M_betaFct                ( 0 ),
    M_count                  ( 0 ),
    M_verbose                ( M_me == 0),
    M_updated                ( false ),
    M_reusePrec              ( false ),
    M_maxIterSolver          ( -1 ),
    M_recomputeMatrix        ( false ),
    M_u                      ( dim_u() ),
    M_p                      ( dim_p() )
{
    M_stab = (&M_uFESpace.refFE() == &M_pFESpace.refFE());
}



template<typename Mesh, typename SolverType>
ADRSolver<Mesh, SolverType>::
~ADRSolver()
{
    delete M_matrMass;
    delete M_matrStokes;
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::setUp( const GetPot& dataFile )
{
    M_steady      = dataFile( "fluid/miscellaneous/steady", 1 );
    M_gammaBeta   = dataFile( "fluid/ipstab/gammaBeta", 0. );
    M_gammaDiv    = dataFile( "fluid/ipstab/gammaDiv", 0. );
    M_gammaPress  = dataFile( "fluid/ipstab/gammaPress", 0. );
    M_divBetaUv   = dataFile( "fluid/space_discretization/div_beta_u_v", 0);
    M_diagonalize = dataFile( "fluid/space_discretization/diagonalize", 1.);

    M_linearSolver.setDataFromGetPot( dataFile, "fluid" );

    M_ipStab.setGammaBeta (M_gammaBeta);
    M_ipStab.setGammaDiv  (M_gammaDiv);
    M_ipStab.setGammaPress(M_gammaPress);

    M_maxIterSolver = dataFile( "fluid/aztec/max_iter", -1);
}

template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::buildSystem()
{

    delete M_matrMass;
    delete M_matrStokes;

    M_matrMass   = new matrix_type(M_localMap);
    M_matrStokes = new matrix_type(M_localMap);

    M_matrMass->zeros();
    M_matrStokes->zeros();

//     if (M_verbose) std::cout << "F-  Velocity unknowns: " << nDimensions*dim_u()     << std::endl;
//     if (M_verbose) std::cout << "F-  Pressure unknowns: " << dim_p()     << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    if (M_verbose) std::cout << "  f-  Computing constant matrices ...        ";

    Chrono chrono;

    Chrono chronoDer;
    Chrono chronoStiff;
    Chrono chronoMass;
    Chrono chronoGrad;

    Chrono chrono1;
    Chrono chrono2;
    Chrono chrono3;
    Chrono chrono4;
    Chrono chrono5;

    // Number of velocity components
    UInt nbCompU = this->u().nbcomp();

    Real bdfCoeff = M_bdf.bdf_u().coeff_der( 0 ) / M_dataType.dataTime()->getTimeStep();
    // Elementary computation and matrix assembling
    // Loop on elements

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();
    UInt pressTotalDof = M_pFESpace.dof().numTotalDof();


    chrono.start();

    for ( UInt iVol = 1; iVol <= M_dataType.dataMesh()->mesh().numElements(); iVol++ )
    {
        chronoDer.start();
        M_pFESpace.fe().update( M_dataType.dataMesh()->mesh().element( iVol ) ); // just to provide the id number in the assem_mat_mixed
//        M_pFESpace.fe().updateFirstDeriv( M_dataType.dataMesh()->mesh().element( iVol ) ); // just to provide the id number in the assem_mat_mixed
        M_uFESpace.fe().updateFirstDeriv( M_dataType.dataMesh()->mesh().element( iVol ) );

        M_elmatStiff.zero();
        M_elmatBdfMass.zero();
        M_elmatMass.zero();
        M_elmatP.zero();
        M_elmatDiv.zero();
        M_elmatGrad.zero();
        chronoDer.stop();


        // stiffness strain
        chronoStiff.start();
        stiff_strain( 2.0*M_dataType.viscosity(), M_elmatStiff, M_uFESpace.fe() );
        //stiff( M_dataType.viscosity(), M_elmatStiff,  M_uFESpace.fe(), 0, 0, nDimensions );
        //stiff_div( 0.5*M_uFESpace.fe().diameter(), M_elmatStiff, M_uFESpace.fe() );
        chronoStiff.stop();

        // mass
        if ( !M_steady )
        {
            chronoMass.start();
            mass( M_dataType.density(), M_elmatMass, M_uFESpace.fe(), 0, 0, nDimensions );
            chronoMass.stop();
        }

        for ( UInt iComp = 0; iComp < nbCompU; iComp++ )
        {
            for ( UInt jComp = 0; jComp < nbCompU; jComp++ )
            {
                // stiffness

                chrono1.start();
                assembleMatrix( *M_matrStokes,
                                M_elmatStiff,
                                M_uFESpace.fe(),
                                M_uFESpace.fe(),
                                M_uFESpace.dof(),
                                M_uFESpace.dof(),
                                iComp, jComp,
                                iComp*velTotalDof, jComp*velTotalDof);
                chrono1.stop();

            }

            if ( !M_steady )
            {
                chrono2.start();
                assembleMatrix( *M_matrMass,
                                M_elmatMass,
                                M_uFESpace.fe(),
                                M_uFESpace.fe(),
                                M_uFESpace.dof(),
                                M_uFESpace.dof(),
                                iComp, iComp,
                                iComp*velTotalDof, iComp*velTotalDof);
                chrono2.stop();
            }

            // div
            chronoGrad.start();
            grad( iComp,  1.0, M_elmatGrad, M_uFESpace.fe(), M_pFESpace.fe(), iComp,     0 );
            chronoGrad.stop();

            chrono3.start();
            assembleMatrix( *M_matrStokes,
                            M_elmatGrad,
                            M_uFESpace.fe(),
                            M_pFESpace.fe(),
                            M_uFESpace.dof(),
                            M_pFESpace.dof(),
                            iComp, 0,
                            iComp*velTotalDof, velTotalDof*nbCompU
                            );
            chrono3.stop();

            chrono4.start();
            assembleTransposeMatrix( *M_matrStokes,
                                     -1.,
                                     M_elmatGrad,
                                     M_pFESpace.fe(),
                                     M_uFESpace.fe(),
                                     M_pFESpace.dof(),
                                     M_uFESpace.dof(),
                                     0 , iComp,
                                     velTotalDof*nbCompU, iComp*velTotalDof
                                     );
            chrono4.stop();
        }
    }

    if (M_stab)
        {
            chrono5.start();
            if (M_verbose)
                std::cout << "  f-  Updating convective ip terms...          "
                          << std::flush;
            vector_type betaVec( *M_localMap.getRepeatedEpetra_Map() );
            M_ipStab.applyPressure( *M_matrStokes, betaVec, M_verbose );

            MPI_Barrier(MPI_COMM_WORLD);

            chrono5.stop();

            if (M_verbose)
                std::cout << "      done in " << chrono5.diff() <<" s." << std::endl << std::flush;
        }


    for (UInt ii = nDimensions*dim_u(); ii < nDimensions*dim_u() + dim_p(); ++ii)
        M_matrStokes->set_mat_inc( ii ,ii, 0. );


    MPI_Barrier(MPI_COMM_WORLD);

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s." << std::endl << std::flush;


    if (M_verbose) std::cout << "  f-  Finalizing the matrices     ...        " << std::flush;
    chrono.start();

    M_matrStokes->getEpetraMatrix().GlobalAssemble();
    M_matrMass->getEpetraMatrix().GlobalAssemble();

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s." << std::endl << std::flush;

}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
updateSystem()
{

    Chrono chrono;

    if (M_recomputeMatrix)
        buildSystem();

    if (M_verbose)
         std::cout << "  f-  Copying constant matrix...               "
                   << std::flush;

    chrono.start();

    Real bdfCoeff = M_bdf.bdf_u().coeff_der( 0 ) / M_dataType.dataTime()->getTimeStep();

    M_matrFull  = new matrix_type(M_localMap);

    *M_matrFull += *M_matrStokes;
    *M_matrFull += *M_matrMass*bdfCoeff;

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s."
                             << std::endl;

    if (M_verbose)
         std::cout << "  f-  Updating the convective terms...         "
                   << std::flush;

    chrono.start();
    UInt nbCompU       = nDimensions;
    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();

    // velocity vector for linearization of convective term

    vector_type betaVec( *M_localMap.getRepeatedEpetra_Map() );

    if ( M_betaFct )
        {
//            this->uInterpolate( *M_betaFct, betaVec, time );
        }
    else
    {
        if ( M_steady )
        {
            betaVec = M_sol; // last iteration
            //betaVec = this->u(); // last iteration
        }
        else
        {
            betaVec = M_bdf.bdf_u().extrap(); // bdf extrapolation
        }
    }


    for ( UInt iVol = 1; iVol<= M_dataType.dataMesh()->mesh().numElements(); ++iVol )
    {

        M_pFESpace.fe().updateFirstDeriv( M_dataType.dataMesh()->mesh().element( iVol ) ); // just to provide the id number in the assem_mat_mixed
        M_uFESpace.fe().updateFirstDeriv( M_dataType.dataMesh()->mesh().element( iVol ) ); //as updateFirstDer

        M_elmatStiff.zero();

        UInt eleID = M_uFESpace.fe().currentLocalId();
        // Non linear term, Semi-implicit approach
        // M_elvec contains the velocity values in the nodes
        for ( UInt iNode = 0 ; iNode < M_uFESpace.fe().nbNode ; iNode++ )
        {
            UInt  iloc = M_uFESpace.fe().patternFirst( iNode );
            for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
            {
                UInt ig = M_uFESpace.dof().localToGlobal( eleID, iloc + 1 ) - 1 + iComp*dim_u();
                M_elvec.vec()[ iloc + iComp*M_uFESpace.fe().nbNode ] = M_dataType.density() * betaVec[ig + 1];
            }
        }

        // Stabilising term: div u^n u v
        if ( M_divBetaUv )
        {
            mass_divw( 0.5*M_dataType.density(), M_elvec, M_elmatStiff, M_uFESpace.fe(), 0, 0, nbCompU );
        }

        // loop on components
        for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
        {
            // compute local convective term and assembling
            grad( 0, M_elvec, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );
            grad( 1, M_elvec, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );
            grad( 2, M_elvec, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );

            assembleMatrix( *M_matrFull,
                            M_elmatStiff,
                            M_uFESpace.fe(),
                            M_uFESpace.fe(),
                            M_uFESpace.dof(),
                            M_uFESpace.dof(),
                            iComp, iComp,
                            iComp*velTotalDof, iComp*velTotalDof
                            );

        }

    }

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s."
                             << std::endl;

    if (M_stab)
    {
        M_ipStab.applyVelocity( *M_matrFull, betaVec, M_verbose );
    }

    M_updated = true;
}



template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
timeAdvance( source_type const& source, Real const& time )
{
//     M_dataType.setTime(time);

    // Number of velocity components
    UInt nbCompU = this->u().nbcomp();

    UInt velTotalDof = M_uFESpace.dof().numTotalDof();

    if (M_verbose)
            std::cout << "  f-  Updating mass term on right hand side... "
                  << std::flush;

    Chrono chrono;

    chrono.start();

    // Right hand side for the velocity at time
    M_rhsNoBC *= 0.;

    // loop on volumes: assembling source term
    for ( UInt iVol = 1; iVol<= M_dataType.dataMesh()->mesh().numElements(); ++iVol )
    {
        M_elvec.zero();
        M_uFESpace.fe().updateJacQuadPt( M_dataType.dataMesh()->mesh().element( iVol ) );

        for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
        {
            // compute local vector
            compute_vec( source, M_elvec, M_uFESpace.fe(), M_dataType.dataTime()->getTime(), iComp );

            // assemble local vector into global one
            assembleVector( M_rhsNoBC, M_elvec,
                            M_uFESpace.fe(), M_uFESpace.dof(),
                            iComp, iComp*velTotalDof );
        }
    }

    if ( !M_steady )
    {
         M_rhsNoBC += *M_matrMass * M_bdf.bdf_u().time_der( M_dataType.dataTime()->getTimeStep() );
    }

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s."  << std::endl;


    M_updated = false;

} // timeAdvance()



template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::iterate( const Real& time )
{

    Chrono chrono;


    // updating the matrix ...

//     if (M_verbose)
//         std::cout << "  f-  Updating the matrix ...    ";

    chrono.start();

    updateSystem();

    chrono.stop();

//     if (M_verbose) std::cout << "done in " << chrono.diff() << " s."
//                              << std::endl;

    // matrix and vector assembling communication
    if (M_verbose)
        std::cout << "  f-  Finalizing the matrix and vectors ...    ";

    chrono.start();


    M_matrFull->getEpetraMatrix().GlobalAssemble();
    M_rhsFull.getEpetraVector().GlobalAssemble();

    chrono.stop();

    if (!M_me)
        std::cout << "done in " << chrono.diff() << " s."
                  << std::endl;

    // boundary conditions update
    if (M_verbose) std::cout << "  f-  Applying boundary conditions...          "
              << std::flush;

    chrono.start();

    applyBoundaryConditions(*M_matrFull);

    chrono.stop();

    if (M_verbose) std::cout << "done in " << chrono.diff() << " s." << std::endl;

    // solving the system

    M_sol = M_bdf.bdf_u().extrap();

    solveSystem();

    M_bdf.bdf_u().shift_right( M_sol );

    delete M_matrFull;
    M_matrFull = 0;

} // iterate()


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::initialize( const Function& x0,
                                                       Real            t0,
                                                       Real            dt )
{
    ID nbComp = this->u().nbcomp(); // Number of components of the velocity

    M_bdf.bdf_u().initialize_unk( x0,
                                  M_dataType.dataMesh()->mesh(),
                                  this->refFEu(),
                                  M_uFESpace.fe(),
                                  M_uFESpace.dof(),
                                  M_localMap,
                                  t0,
                                  dt,
                                  nbComp + 1 );

    M_sol = *( M_bdf.bdf_u().unk().begin() );

    vector_type vel(M_sol, 0);

    if (M_me == 0)
    {
            for ( UInt iDof = 0; iDof<nDimensions*dim_u(); ++iDof )
            {
                this->u()[ iDof ] = vel[ iDof + 1 ];
            }

            for ( UInt iDof = 0; iDof<dim_p(); ++iDof )
            {
                this->p()[ iDof ] = vel[ iDof + nDimensions*dim_u() + 1 ];
            }

    }


} // initialize()

/*! Initialize when initial values for the velocity and the pressure are read
  from file
  @author Martin Prosi
*/

template<typename Mesh, typename SolverType>
void
ADRSolver<Mesh, SolverType>::initialize( const std::string & vname )
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

    for ( UInt iDof = 0; iDof < nDimensions*dim_u(); ++iDof )
    {
        this->u()[ iDof ] = M_sol[ iDof ];
    }
    for ( UInt iDof = 0; iDof<dim_u(); ++iDof )
    {
        this->p()[ iDof ] = M_sol[ iDof+nDimensions*dim_u() ];
    }

    M_bdf.bdf_u().initialize_unk( M_sol );
    //bdf().bdf_p().initialize_unk( _p() );

    M_bdf.bdf_u().showMe();
    //bdf().bdf_p().showMe();

}

template<typename Mesh, typename SolverType>
void
ADRSolver<Mesh, SolverType>::initializeNull()
{
    vector_type null(M_localMap);
    M_bdf.bdf_u().initialize_unk( M_sol );
}


template<typename Mesh, typename SolverType>
void
ADRSolver<Mesh, SolverType>::initializeStokes( source_type const& source,
                                                        Real t0 )
{

    bool isSteady = M_steady;

    M_steady      = true;
    timeAdvance( source, t0 );
    M_steady      = isSteady;

    int NumProc   = M_comm->NumProc();

    Chrono chrono;

    // velocity vector for linearization of convective term
    vector_type betaVec( *M_localMap.getRepeatedEpetra_Map() );

//    M_dataType.setTime(t0);


    if (M_verbose) std::cout << "  f-  Copying Stokes matrix...                 "
              << std::flush;
    chrono.start();

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();


//    matrix_type matrFull( M_localMap );
    M_matrFull = new matrix_type( M_localMap );
    *M_matrFull += *M_matrStokes;

    chrono.stop();

    if (M_verbose) std::cout << "done in " << chrono.diff() << " s."
                           << std::endl;


#if 1
    if (M_stab)
    {
        if (M_verbose) std::cout << "  f-  Adding ip terms...                       "
                               << std::flush;
        chrono.start();
//        initStab.applyPressure( matrFull, betaVec );
        M_ipStab.applyVelocity( *M_matrFull, betaVec );
        chrono.stop();

        if (M_verbose) std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
#endif

    if (M_verbose)
        std::cout << "  f-  Finalizing the matrix and vectors ...    ";

    chrono.start();

    M_matrFull->getEpetraMatrix().GlobalAssemble();
//    M_rhsFull.getEpetraVector().GlobalAssemble();

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s."
                  << std::endl;


    if (M_verbose) std::cout << "  f-  Applying boundary conditions   ...       "
                  << std::flush;

    chrono.start();
    applyBoundaryConditions(*M_matrFull);
    chrono.stop();

    MPI_Barrier(MPI_COMM_WORLD);

    if (M_verbose)
    std::cout << "done in " << chrono.diff() << " s."
              << std::endl;

    solveSystem();

    M_bdf.bdf_u().initialize_unk( M_sol );

    delete M_matrFull;
    M_matrFull = 0;

    M_linearSolver.precReset();
} // initializeStokes


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::solveSystem()
{

}



template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::removeMean( vector_type& x, UInt comp )
{
    Real sum1 = 0.;
    Real sum0 = 0.;

    for ( UInt iVol = 1; iVol <= M_dataType.dataMesh()->mesh().numElements(); iVol++ )
    {
        M_pFESpace.fe().updateFirstDeriv( M_dataType.dataMesh()->mesh().element( iVol ) );
        sum1 += elem_integral( x, M_pFESpace.fe(), M_pFESpace.dof(), comp );
        sum0 += M_pFESpace.fe().measure();
    }
    Real mean = sum1/sum0;

    std::cout << "   Mean pressure value : " << mean << std::endl;
    for ( UInt iNode = (comp - 1)*dim_u(); iNode < comp*dim_u(); ++iNode )
    {
        x[ iNode + 1 ] -= mean;
    }

} // removeMean()

template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::applyBoundaryConditions(matrix_type &_csr)
{
    M_rhsFull = M_rhsNoBC;

    // BC manage for the velocity
    if ( !this->bcHandler().bdUpdateDone() )
        this->bcHandler().bdUpdate( M_dataType.dataMesh()->mesh(), M_uFESpace.feBd(), M_uFESpace.dof() );

    bcManage( _csr, M_rhsFull, M_dataType.dataMesh()->mesh(), M_uFESpace.dof(), this->bcHandler(), M_uFESpace.feBd(), 1.,
               M_dataType.dataTime()->getTime() );

    if ( this->bcHandler().hasOnlyEssential() && M_diagonalize )
    {
        _csr.diagonalize( nDimensions*dim_u(),
                          M_diagonalize,
                          M_rhsFull,
                          0.);
    }
} // applyBoundaryCondition




// Postprocessing
template <typename Mesh, typename SolverType>
void
ADRSolver<Mesh, SolverType>::postProcess()
{
    std::ostringstream index;
    std::string name;


    if ( fmod( float( M_count ), float( M_dataType.verbose() ) ) == 0.0 )
    {
        std::cout << "  F-  Post-processing \n";
        index << std::setfill('0') << std::setw(3);
        index << ( M_count / M_dataType.verbose() );
        name = index.str();

        // postprocess data file for medit
        wr_medit_ascii_scalar( "press." + name + ".bb", this->p().giveVec(),
                               this->p().size() );
        wr_medit_ascii_scalar( "vel_x." + name + ".bb", this->u().giveVec(),
                               M_dataType.dataMesh()->mesh().numGlobalVertices() );
        wr_medit_ascii_scalar( "vel_y." + name + ".bb", this->u().giveVec() + this->dim_u(),
                               M_dataType.dataMesh()->mesh().numGlobalVertices() );
        wr_medit_ascii_scalar( "vel_z." + name + ".bb", this->u().giveVec()+2*this->dim_u(),
                               M_dataType.dataMesh()->mesh().numGlobalVertices() );
        system( ( "ln -s -f " + M_dataType.dataMesh()->meshDir() + M_dataType.dataMesh()->meshFile() +
                  " press." + name + ".mesh" ).data() );
        system( ( "ln -s -f " + M_dataType.dataMesh()->meshDir() + M_dataType.dataMesh()->meshFile() +
                  " vel_x." + name + ".mesh" ).data() );
        system( ( "ln -s -f " + M_dataType.dataMesh()->meshDir() + M_dataType.dataMesh()->meshFile() +
                  " vel_y." + name + ".mesh" ).data() );
        system( ( "ln -s -f " + M_dataType.dataMesh()->meshDir() + M_dataType.dataMesh()->meshFile() +
                  " vel_z." + name + ".mesh" ).data() );
    }

    M_count++;
}



template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::echo(std::string message)
{
    if (!M_me) std::cout << message << std::flush;
}

} // namespace LifeV


#endif //_ADRSOLVER_H_
