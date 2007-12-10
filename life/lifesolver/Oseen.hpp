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
  \file Oseen.hpp
  \author G. Fourestey
  \date 2/2007

  \brief This file contains a Oseen equation solver class
*/
#ifndef _OSEEN_H_
#define _OSEEN_H_

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifearray/pattern.hpp>
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

#include "life/lifefem/FESpace.hpp"

namespace LifeV
{
/*!
  \class Oseen

  This class implements a NavierStokes solver.
  The resulting linear systems are solved by GMRES on the full
  matrix ( u and p coupled ).

*/



template< typename Mesh,
          typename SolverType = LifeV::Epetra::SolverTrilinos >
class Oseen
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
    typedef boost::shared_ptr<matrix_type>        matrix_ptrtype;
    typedef typename SolverType::vector_type      vector_type;

    typedef typename SolverType::prec_raw_type    prec_raw_type;
    typedef typename SolverType::prec_type        prec_type;


    //! Constructor
    /*!
      \param dataType
      \param localMap
      \param velocity FE space
      \param pressure FE space
      \param bcHandler boundary conditions for the velocity
    */
    Oseen( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& uFESpace,
           FESpace<Mesh, EpetraMap>& pFESpace,
           BCHandler&                bcHandler,
           Epetra_Comm&              comm );

    /*!
      \param dataType
      \param localMap
      \param velocity FE space
      \param pressure FE space
    */
    Oseen( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& uFESpace,
           FESpace<Mesh, EpetraMap>& pFESpace,
           Epetra_Comm&              comm );

    //! virtual destructor

    virtual ~Oseen();

    //! Update convective term, bc treatment and solve the linearized ns system

    virtual void iterate( bchandler_raw_type& bch );

    virtual void setUp        ( const GetPot& dataFile );

    virtual void buildSystem();

    virtual void updateRHS(vector_type const& sourceVec)
    {
        M_rhsNoBC = sourceVec;
        M_rhsNoBC.GlobalAssemble();
    }

    virtual void updateSystem(double       alpha,
                              vector_type& betaVec,
                              vector_type& sourceVec
                              );

//     void updateLinearSystem(double       alpha,
//                             vector_type& betaVec,
//                             vector_type& sourceVec
//                             );

    void initialize( const Function&, const Function&  );
    void initialize( const vector_type&, const vector_type& );

    //! returns the local solution vector

    const vector_type& solution() const {return M_sol;}

    //! returns the local residual vector

    const vector_type& residual() const {return M_residual;}

    //! reduce the local solution solution in global vectors

    void reduceSolution( Vector& u,
                         Vector& p );

    void reduceResidual( Vector& res);


    FESpace<Mesh, EpetraMap>& velFESpace()   {return M_uFESpace;}
    FESpace<Mesh, EpetraMap>& pressFESpace() {return M_pFESpace;}


    //! Bounday Conditions
    const bool BCset() const {return M_setBC;}
    //! set the fluid BCs
    void setBC(BCHandler &BCh_u)
        {
            M_BCh_fluid = &BCh_u; M_setBC = true;
        }

//     BCHandler& bcHandler()
//         {
//             return *M_BCh_fluid;
//         }

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
    void postProcess(bool _writeMesh = false);

    void resetPrec() {M_resetPrec = true;}

    Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.getRepeatedEpetra_Map(); }

    EpetraMap const& getMap() const { return M_localMap; }

    const Epetra_Comm& comm() const {return *M_comm;}

    void recomputeMatrix(bool const recomp){M_recomputeMatrix = recomp;}

    matrix_type& matrMass()
        {
            return *M_matrMass;
        }


protected:


    void solveSystem            (  matrix_ptrtype matrFull,
                                   vector_type&   rhsFull );

    void applyBoundaryConditions(  matrix_type&        matrix,
                                   vector_type&        rhs,
                                   bchandler_raw_type& BCh);

    void echo(std::string message);

    //! removes mean of component comp of vector x
    void removeMean( vector_type& x, UInt comp = 1 );

    //private members

    //! data for NS solvers
    const data_type&               M_data;

    // FE spaces
    FESpace<Mesh, EpetraMap>&      M_uFESpace;
    FESpace<Mesh, EpetraMap>&      M_pFESpace;

    //! MPI communicator
    Epetra_Comm*                   M_comm;
    int                            M_me;

    //! fluid BC
    BCHandler*                     M_BCh_fluid;
    bool                           M_setBC;

    EpetraMap                      M_localMap;

    //! mass matrix

    matrix_ptrtype                   M_matrMass;

    //! Stokes matrix: mu*stiff
    matrix_ptrtype                   M_matrStokes;

    //! matrix to be solved
//    matrix_ptrtype                 M_matrFull;

    //! matrix without boundary conditions

    matrix_ptrtype                 M_matrNoBC;

    //! source term for NS
    source_type                    M_source;


    //! Right hand side for the velocity
    vector_type                    M_rhsNoBC;

    //! Right hand side global
    vector_type                    M_rhsFull;

    //! Global solution _u and _p
    vector_type                    M_sol;

    //! residual

    vector_type                    M_residual;

    SolverType                     M_linearSolver;

    prec_type                      M_prec;

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
    int                            M_maxIterForReuse;
    bool                           M_resetPrec;

    //! interger storing the max number of solver iteration with prec recomputing

    int                            M_maxIterSolver;

    //!

    bool                           M_recomputeMatrix;

private:

    //! Elementary matrices and vectors
    ElemMat                        M_elmatStiff;      // velocity Stokes
    ElemMat                        M_elmatMass;       // velocity mass
    ElemMat                        M_elmatP;          // (p,q) bloc for preconditioners
    ElemMat                        M_elmatDiv;
    ElemMat                        M_elmatGrad;
    ElemVec                        M_elvec;           // Elementary right hand side

    UInt dim_u() const           { return M_uFESpace.dim(); }
    UInt dim_p() const           { return M_pFESpace.dim(); }

}; // class Oseen



//
// IMPLEMENTATION
//


template<typename Mesh, typename SolverType>
Oseen<Mesh, SolverType>::
Oseen( const data_type&          dataType,
       FESpace<Mesh, EpetraMap>& uFESpace,
       FESpace<Mesh, EpetraMap>& pFESpace,
       BCHandler&                BCh_u,
       Epetra_Comm&              comm ):
    M_data                   ( dataType ),
    M_uFESpace               ( uFESpace ),
    M_pFESpace               ( pFESpace ),
    M_BCh_fluid              ( &BCh_u ),
    M_setBC                  ( true ),
    M_comm                   ( &comm ),
    M_me                     ( M_comm->MyPID() ),
    M_linearSolver           ( ),
    M_prec                   ( new prec_raw_type() ),
    M_localMap               ( M_uFESpace.map()  + M_pFESpace.map() ),
    M_matrMass               ( ),
    M_matrStokes             ( ),
//    M_matrFull               ( ),
    M_matrNoBC               ( ),
    M_elmatStiff             ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatMass              ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatP                 ( M_pFESpace.fe().nbNode, 1, 1 ),
    M_elmatDiv               ( M_pFESpace.fe().nbNode, 1, 0, M_uFESpace.fe().nbNode, 0, nDimensions ),
    M_elmatGrad              ( M_uFESpace.fe().nbNode, nDimensions, 0, M_pFESpace.fe().nbNode, 0, 1 ),
    M_elvec                  ( M_uFESpace.fe().nbNode, nDimensions ),
    M_rhsNoBC                ( M_localMap ),
    M_rhsFull                ( M_localMap ),
    M_sol                    ( M_localMap ),
    M_residual               ( M_localMap ),
    M_stab                   ( false ),
    M_ipStab                 ( M_uFESpace.mesh(),
                               M_uFESpace.dof(), M_uFESpace.refFE(),
                               M_uFESpace.feBd(), M_uFESpace.qr(),
                               0., 0., 0.,
                               M_data.viscosity() ),
    M_betaFct                ( 0 ),
    M_count                  ( 0 ),
    M_verbose                ( M_me == 0),
    M_updated                ( false ),
    M_reusePrec              ( true ),
    M_resetPrec              ( true ),
    M_maxIterSolver          ( -1 ),
    M_recomputeMatrix        ( false )
{
    M_stab = (&M_uFESpace.refFE() == &M_pFESpace.refFE());
}



template<typename Mesh, typename SolverType>
Oseen<Mesh, SolverType>::
Oseen( const data_type&          dataType,
       FESpace<Mesh, EpetraMap>& uFESpace,
       FESpace<Mesh, EpetraMap>& pFESpace,
       Epetra_Comm&              comm ):
    M_data               ( dataType ),
    M_uFESpace               ( uFESpace ),
    M_pFESpace               ( pFESpace ),
    M_setBC                  ( false ),
    M_comm                   ( &comm ),
    M_me                     ( M_comm->MyPID() ),
    M_linearSolver           ( ),
    M_prec                   ( new prec_raw_type() ),
    M_localMap               ( M_uFESpace.map() + M_pFESpace.map() ),
    M_matrMass               ( ),
    M_matrStokes             ( ),
//    M_matrFull               ( ),
    M_matrNoBC               ( ),
    M_elmatStiff             ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatMass              ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatP                 ( M_pFESpace.fe().nbNode, 1, 1 ),
    M_elmatDiv               ( M_pFESpace.fe().nbNode, 1, 0, M_uFESpace.fe().nbNode, 0, nDimensions ),
    M_elmatGrad              ( M_uFESpace.fe().nbNode, nDimensions, 0, M_pFESpace.fe().nbNode, 0, 1 ),
    M_elvec                  ( M_uFESpace.fe().nbNode, nDimensions ),
    M_rhsNoBC                ( M_localMap ),
    M_residual               ( M_localMap ),
    M_rhsFull                ( M_localMap ),
    M_sol                    ( M_localMap ),
    M_stab                   ( false ),
    M_ipStab                 ( M_data.mesh(),
                               M_uFESpace.dof(), M_uFESpace.refFE(),
                               M_uFESpace.feBd(), M_uFESpace.qr(),
                               0., 0., 0.,
                               M_data.viscosity() ),
    M_betaFct                ( 0 ),
    M_count                  ( 0 ),
    M_verbose                ( M_me == 0),
    M_updated                ( false ),
    M_reusePrec              ( true ),
    M_maxIterForReuse        ( -1 ),
    M_resetPrec              ( true ),
    M_maxIterSolver          ( -1 ),
    M_recomputeMatrix        ( false )
{
    M_stab = (&M_uFESpace.refFE() == &M_pFESpace.refFE());
}



template<typename Mesh, typename SolverType>
Oseen<Mesh, SolverType>::
~Oseen()
{

}


template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::setUp( const GetPot& dataFile )
{
    M_steady      = dataFile( "fluid/miscellaneous/steady",        1  );
    M_gammaBeta   = dataFile( "fluid/ipstab/gammaBeta",            0. );
    M_gammaDiv    = dataFile( "fluid/ipstab/gammaDiv",             0. );
    M_gammaPress  = dataFile( "fluid/ipstab/gammaPress",           0. );
    M_divBetaUv   = dataFile( "fluid/discretization/div_beta_u_v", 0  );
    M_diagonalize = dataFile( "fluid/discretization/diagonalize",  1. );


    M_linearSolver.setDataFromGetPot( dataFile, "fluid/solver" );

    M_ipStab.setGammaBeta (M_gammaBeta);
    M_ipStab.setGammaDiv  (M_gammaDiv);
    M_ipStab.setGammaPress(M_gammaPress);

    M_maxIterSolver   = dataFile( "fluid/solver/max_iter", -1);
    M_reusePrec       = dataFile( "fluid/prec/reuse", true);
    M_maxIterForReuse = dataFile( "fluid/prec/max_iter_reuse", M_maxIterSolver*8/10);

    M_prec->setDataFromGetPot( dataFile, "fluid/prec" );
}

template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::buildSystem()
{


    M_matrMass.reset  ( new matrix_type(M_localMap) );
    M_matrStokes.reset( new matrix_type(M_localMap) );

//    M_comm->Barrier();

    if (M_verbose) std::cout << "  f-  Computing constant matrices ...        ";

    Chrono chrono;

    Chrono chronoDer;
    Chrono chronoStiff;
    Chrono chronoMass;
    Chrono chronoGrad;

    Chrono chronoStiffAssemble;
    Chrono chronoMassAssemble;
    Chrono chronoGradAssemble;
    Chrono chronoDivAssemble;
    Chrono chronoStab;
    Chrono chronoZero;

    // Number of velocity components
    UInt nbCompU = nDimensions;

    // Elementary computation and matrix assembling
    // Loop on elements

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();
//    UInt pressTotalDof = M_pFESpace.dof().numTotalDof();


    chrono.start();

    for ( UInt iVol = 1; iVol <= M_uFESpace.mesh()->numVolumes(); iVol++ )
    {
        chronoDer.start();
        M_pFESpace.fe().update( M_uFESpace.mesh()->volumeList( iVol ) ); // just to provide the id number in the assem_mat_mixed
//        M_pFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) ); // just to provide the id number in the assem_mat_mixed
        M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) );

        chronoDer.stop();

        chronoZero.start();
        M_elmatStiff.zero();
        M_elmatMass.zero();
        M_elmatP.zero();
        M_elmatDiv.zero();
        M_elmatGrad.zero();
        chronoZero.stop();


        // stiffness strain
        chronoStiff.start();
        //stiff_strain( 2.0*M_data.viscosity(), M_elmatStiff, M_uFESpace.fe() );
        stiff( M_data.viscosity(), M_elmatStiff,  M_uFESpace.fe(), 0, 0, nDimensions );
        //stiff_div( 0.5*M_uFESpace.fe().diameter(), M_elmatStiff, M_uFESpace.fe() );
        chronoStiff.stop();

        // mass
        if ( !M_steady )
        {
            chronoMass.start();
            mass( M_data.density(), M_elmatMass, M_uFESpace.fe(), 0, 0, nDimensions );
            chronoMass.stop();
        }

        for ( UInt iComp = 0; iComp < nbCompU; iComp++ )
        {
//             for ( UInt jComp = 0; jComp < nbCompU; jComp++ )
//             {
                // stiffness

                chronoStiffAssemble.start();
                assembleMatrix( *M_matrStokes,
                                M_elmatStiff,
                                M_uFESpace.fe(),
                                M_uFESpace.fe(),
                                M_uFESpace.dof(),
                                M_uFESpace.dof(),
                                iComp, iComp,
                                iComp*velTotalDof, iComp*velTotalDof);
                chronoStiffAssemble.stop();

//             }

            if ( !M_steady )
            {
                chronoMassAssemble.start();
                assembleMatrix( *M_matrMass,
                                M_elmatMass,
                                M_uFESpace.fe(),
                                M_uFESpace.fe(),
                                M_uFESpace.dof(),
                                M_uFESpace.dof(),
                                iComp, iComp,
                                iComp*velTotalDof, iComp*velTotalDof);
                chronoMassAssemble.stop();
            }

            // div
            chronoGrad.start();
            grad( iComp,  1.0, M_elmatGrad, M_uFESpace.fe(), M_pFESpace.fe(), iComp,     0 );
            chronoGrad.stop();

            chronoGradAssemble.start();
            assembleMatrix( *M_matrStokes,
                            M_elmatGrad,
                            M_uFESpace.fe(),
                            M_pFESpace.fe(),
                            M_uFESpace.dof(),
                            M_pFESpace.dof(),
                            iComp, 0,
                            iComp*velTotalDof, velTotalDof*nbCompU
                            );
            chronoGradAssemble.stop();

            chronoDivAssemble.start();
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
            chronoDivAssemble.stop();
        }
    }

//     if (M_stab)
//         {
//             chronoStab.start();
//             if (M_verbose)
//                 std::cout << "  f-  Updating convective ip terms... "
//                           << std::flush;
//             vector_type betaVec( *M_localMap.getRepeatedEpetra_Map() );
//             if (M_verbose)
//                 std::cout << "(repeated vector created in " << chronoStab.diff() <<" s.)" << std::flush;
//             M_ipStab.applyPressure( *M_matrStokes, betaVec, M_verbose );
//             M_ipStab.applyVelocity( *M_matrStokes, betaVec, M_verbose );

//             //MPI_Barrier(MPI_COMM_WORLD);

//             chronoStab.stop();

//             if (M_verbose)
//                 std::cout << "      done in " << chronoStab.diff() <<" s.\n" << std::flush;
//         }


    for (UInt ii = nDimensions*dim_u(); ii < nDimensions*dim_u() + dim_p(); ++ii)
        M_matrStokes->set_mat_inc( ii ,ii, 0. );


    M_comm->Barrier();

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;


    if (M_verbose) std::cout << "  f-  Finalizing the matrices     ...        " << std::flush;
    chrono.start();

    M_matrStokes->GlobalAssemble();
    M_matrMass->GlobalAssemble();

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s." << std::endl;

    if (false)
        std::cout << "partial times:  \n"
                  << " Der            " << chronoDer.diff_cumul() << " s.\n"
                  << " Stab           " << chronoStab.diff_cumul() << " s.\n"
                  << " Zero           " << chronoZero.diff_cumul() << " s.\n"
                  << " Stiff          " << chronoStiff.diff_cumul() << " s.\n"
                  << " Stiff Assemble " << chronoStiffAssemble.diff_cumul() << " s.\n"
                  << " Mass           " << chronoMass.diff_cumul() << " s.\n"
                  << " Mass Assemble  " << chronoMassAssemble.diff_cumul() << " s.\n"
                  << " Grad           " << chronoGrad.diff_cumul() << " s.\n"
                  << " Grad Assemble  " << chronoGradAssemble.diff_cumul() << " s.\n"
                  << " Div Assemble   " << chronoDivAssemble.diff_cumul() << " s.\n"
                  << std::endl;

}


template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::
initialize( const Function& u0, const Function& p0 )
{
     vector_type u(M_uFESpace.map());
     M_uFESpace.interpolate(u0, u, 0.0);

     vector_type p(M_pFESpace.map());
     M_pFESpace.interpolate(p0, p, 0.0);

     initialize(u, p);
}


template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::
initialize( const vector_type& u0, const vector_type& p0 )
{

    M_sol = u0;

    M_sol.add(p0, nDimensions*M_uFESpace.dof().numTotalDof());

}


template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::
updateSystem(double       alpha,
             vector_type& betaVec,
             vector_type& sourceVec
             )
{

    Chrono chrono;

    if (M_verbose)
            std::cout << "  f-  Updating mass term on right hand side... "
                  << std::flush;

    chrono.start();

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();
//    UInt pressTotalDof = M_pFESpace.dof().numTotalDof();

    // Right hand side for the velocity at time

    updateRHS(sourceVec);

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"  << std::flush;


    M_updated = false;

//

    if (M_recomputeMatrix)
        buildSystem();

    if (M_verbose)
          std::cout << "  f-  Copying the matrices ...                 "
                    << std::flush;

    chrono.start();

    M_matrNoBC.reset(new matrix_type(M_localMap, M_matrStokes->getMeanNumEntries() ));

    *M_matrNoBC += *M_matrStokes;

    if (alpha != 0. )
    {
        *M_matrNoBC += *M_matrMass*alpha;
    }


    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n"
                             << std::flush;


    UInt nbCompU       = nDimensions;

    //! managing the convective term

    double normInf;
    betaVec.NormInf(&normInf);

    if (normInf != 0.)
    {

        if (M_verbose)
            std::cout << "  f-  Updating the convective terms ...        "
                      << std::flush;

        // vector with repeated nodes over the processors

        vector_type betaVecRep(betaVec,*M_localMap.getRepeatedEpetra_Map() );

        chrono.start();

        for ( UInt iVol = 1; iVol<= M_uFESpace.mesh()->numVolumes(); ++iVol )
        {

            M_pFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) ); // just to provide the id number in the assem_mat_mixed
            M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) ); //as updateFirstDer

            M_elmatStiff.zero();

            UInt eleID = M_uFESpace.fe().currentLocalId();
            // Non linear term, Semi-implicit approach
            // M_elvec contains the velocity values in the nodes
            for ( UInt iNode = 0 ; iNode < ( UInt ) M_uFESpace.fe().nbNode ; iNode++ )
            {
                UInt  iloc = M_uFESpace.fe().patternFirst( iNode );
                for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
                {
                    UInt ig = M_uFESpace.dof().localToGlobal( eleID, iloc + 1 ) + iComp*dim_u();
                    M_elvec.vec()[ iloc + iComp*M_uFESpace.fe().nbNode ] = M_data.density() * betaVecRep[ig]; // BASEINDEX + 1
                }
            }

            // Stabilising term: div u^n u v
            if ( M_divBetaUv )
            {
                mass_divw( 0.5*M_data.density(), M_elvec, M_elmatStiff, M_uFESpace.fe(), 0, 0, nbCompU );
            }

            // loop on components
            for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
            {
                // compute local convective term and assembling
                grad( 0, M_elvec, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );
                grad( 1, M_elvec, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );
                grad( 2, M_elvec, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );

                assembleMatrix( *M_matrNoBC,
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
        if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n"
                                 << std::flush;

        if (M_stab)
        {  // if there is a problem of indeces here, you should use the repeated version of beta
            //M_ipStab.applyPressure( *M_matrFull, betaVecRep, M_verbose );
            M_ipStab.applyVelocity( *M_matrNoBC, betaVecRep, M_verbose );
        }

    }

    M_updated = true;

     if (alpha != 0.)
     {
         M_matrNoBC->GlobalAssemble();
     }

}




template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::iterate( bchandler_raw_type& bch )
{

    Chrono chrono;


    // updating the matrix ...

//     if (M_verbose)
//         std::cout << "  f-  Updating the matrix ...    ";

//     chrono.start();

// //    updateSystem);

//     chrono.stop();


    // matrix and vector assembling communication

    if (M_verbose)
        {
            std::cout << "  f-  Finalizing the matrix and vectors ...    ";
        }

    chrono.start();


    M_matrNoBC->GlobalAssemble();

    matrix_ptrtype matrFull( new matrix_type(*M_matrNoBC) );
    vector_type    rhsFull = M_rhsNoBC;

//     matrFull.reset(new matrix_type(*M_matrNoBC));
//     M_rhsFull = M_rhsNoBC;

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    // boundary conditions update
    M_comm->Barrier();
    if (M_verbose) std::cout << "  f-  Applying boundary conditions ...         "
              << std::flush;

    chrono.start();

    applyBoundaryConditions( *matrFull, rhsFull, bch);

    chrono.stop();

    M_comm->Barrier();

    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;
    // solving the system

    solveSystem( matrFull, rhsFull );

    M_residual  = M_rhsNoBC;
    M_residual -= *M_matrNoBC*M_sol;

} // iterate()



template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::solveSystem( matrix_ptrtype  matrFull,
                                           vector_type&    rhsFull )
{
    Chrono chrono;

    if (M_verbose)
        std::cout << "  f-  Setting up the solver ...                ";

    chrono.start();
//    assert (M_matrFull.get() != 0);
    M_linearSolver.setMatrix(*matrFull);
    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    // overlapping schwarz preconditioner

    if ( !M_reusePrec || M_resetPrec || !M_prec->set() )
    {
        chrono.start();

        if (M_verbose)
            std::cout << "  f-  Computing the precond ...                ";

        M_prec->buildPreconditioner(matrFull);

        double condest = M_prec->Condest();

        M_linearSolver.setPreconditioner(M_prec);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "  f-       Estimated condition number = " << condest << "\n" <<  std::flush;
        }


    }
    else
    {
        if (M_verbose)
            std::cout << "  f-  Reusing  precond ...                \n" <<  std::flush;
    }


    chrono.start();

    if (M_verbose)
        std::cout << "  f-  Solving system ...                                ";

    int numIter = M_linearSolver.solve(M_sol, rhsFull);

    if (numIter > M_maxIterSolver)
    {
        chrono.start();

        if (M_verbose)
            std::cout << "  f- Iterative solver failed, recomputing the precond ...                ";

        M_prec->buildPreconditioner(matrFull);

        double condest = M_prec->Condest();

        M_linearSolver.setPreconditioner(M_prec);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "  f-      Estimated condition number = " << condest << "\n" <<  std::flush;
        }

        numIter = M_linearSolver.solve(M_sol, rhsFull);

        if (numIter > M_maxIterSolver && M_verbose)
            std::cout << "  f- ERROR: Iterative solver failed again.\n";

    }

    M_resetPrec = (numIter > M_maxIterForReuse);

    chrono.stop();
    if (M_verbose)
    {
        std::cout << "\ndone in " << chrono.diff()
                  << " s. ( " << numIter << "  iterations. ) \n"
                  << std::flush;
    }

    M_comm->Barrier();

}


template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::reduceSolution( Vector& u,
                                              Vector& p )
{
    vector_type vel(M_sol, 0);

    if (M_verbose)
    {
        for ( UInt iDof = 0; iDof < nDimensions*dim_u(); ++iDof )
        {
            u[ iDof ] = vel[ iDof + 1 ]; // BASEINDEX + 1
        }

        for ( UInt iDof = 0; iDof<dim_p(); ++iDof )
        {
            p[ iDof ] = vel[ iDof + nDimensions*dim_u() + 1 ]; // BASEINDEX + 1
        }
    }

}

template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::reduceResidual( Vector& res )
{
    vector_type vel(M_residual, 0);

    if (M_verbose)
    {
        for ( UInt iDof = 0; iDof < nDimensions*dim_u(); ++iDof )
        {
            res[ iDof ] = vel[ iDof + 1 ]; // BASEINDEX + 1
        }

    }
}


template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::removeMean( vector_type& x, UInt comp )
{
    Real sum1 = 0.;
    Real sum0 = 0.;

    for ( UInt iVol = 1; iVol <= M_uFESpace.mesh()->numVolumes(); iVol++ )
    {
        M_pFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) );
        sum1 += elem_integral( x, M_pFESpace.fe(), M_pFESpace.dof(), comp );
        sum0 += M_pFESpace.fe().measure();
    }
    Real mean = sum1/sum0;

    std::cout << "   Mean pressure value : " << mean << std::endl;
    for ( UInt iNode = (comp - 1)*dim_u(); iNode < comp*dim_u(); ++iNode )
    {
        x[ iNode + 1 ] -= mean; // BASEINDEX + 1
    }

} // removeMean()

template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::applyBoundaryConditions( matrix_type&        matrix,
                                                       vector_type&        rhs,
                                                       bchandler_raw_type& BCh )
{

    // M_rhsFull = M_rhsNoBC;

    // BC manage for the velocity
    if ( !BCh.bdUpdateDone() )
    {
        BCh.bdUpdate( *M_uFESpace.mesh(), M_uFESpace.feBd(), M_uFESpace.dof() );
    }

    vector_type rhsFull(*M_localMap.getRepeatedEpetra_Map());


    rhsFull.Import(M_rhsNoBC, Zero); // ignoring non-local entries, Otherwise they are summed up lately

    bcManage( matrix, rhsFull, *M_uFESpace.mesh(), M_uFESpace.dof(), BCh, M_uFESpace.feBd(), 1.,
              M_data.time() );

    rhs = rhsFull;


    if ( BCh.hasOnlyEssential() && M_diagonalize )
    {
        matrix.diagonalize( nDimensions*dim_u(),
                            M_diagonalize,
                            rhs,
                            0.);
    }

} // applyBoundaryCondition





// Postprocessing
template <typename Mesh, typename SolverType>
void
Oseen<Mesh, SolverType>::postProcess(bool _writeMesh)
{
    std::ostringstream index;
    std::ostringstream indexMe;
    std::string name;
    std::string me;

    M_count++;

    indexMe << M_me;

    switch ( indexMe.str().size() )
    {
        case 1:
            me = "00" + indexMe.str();
            break;
        case 2:
            me = "0" + indexMe.str();
            break;
        case 3:
            me = indexMe.str();
            break;
    }


//     if (_writeMesh || (M_count / M_data.verbose() == 0) )
//         writeMesh  ("partedMesh." + me + ".mesh", M_pFESpace.mesh() );

//     vector_type velAndPressure(M_sol,*M_localMap.getRepeatedEpetra_Map());


//         if ( fmod( float( M_count ), float( M_data.verbose() ) ) == 0.0 )
//         {
    if (M_me == 0)
        std::cout << "  F-  Post-processing " << std::flush;

    index << std::setfill('0') << std::setw(3);
    index << ( M_count / M_data.verbose() );
    name = index.str();

    PhysVectUnknown<Vector> u(nDimensions*dim_u());
    ScalUnknown<Vector>     p(dim_p());

    reduceSolution(u, p);

    if (M_me == 0)
    {
        // postprocess data file for medit
        wr_medit_ascii_scalar( "vel_x." + name + ".bb", u.giveVec(),
                               M_data.mesh()->numGlobalVertices() );
        wr_medit_ascii_scalar( "vel_y." + name + ".bb", u.giveVec() + this->dim_u(),
                               M_data.mesh()->numGlobalVertices() );
        wr_medit_ascii_scalar( "vel_z." + name + ".bb", u.giveVec()+2*this->dim_u(),
                               M_data.mesh()->numGlobalVertices() );
        wr_medit_ascii_scalar( "press." + name + ".bb", p.giveVec(),
                               p.size() );

//  meditSolutionWriter( "vel_x." + name + "." + me + ".bb",
//                              M_uFESpace.mesh(), velAndPressure, M_uFESpace.dof().numTotalDof()*0);
//         meditSolutionWriter( "vel_y." + name + "." + me + ".bb",
//                              M_uFESpace.mesh(), velAndPressure, M_uFESpace.dof().numTotalDof()*1);
//         meditSolutionWriter( "vel_z." + name + "." + me + ".bb",
//                              M_uFESpace.mesh(), velAndPressure, M_uFESpace.dof().numTotalDof()*2);


//         meditSolutionWriter( "press." + name + "." + me + ".bb",
//                              M_pFESpace.mesh(), velAndPressure, M_uFESpace.dof().numTotalDof()*3);


//    cout << ( "should do ln -s -f " + M_data.meshDir() + M_data.meshFile() +
//             " press." + name + "." + me + ".mesh" ).data()  ;

        system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
                  " press." + name + ".mesh" ).data() );
        system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
                  " vel_x." + name + ".mesh" ).data() );
        system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
                  " vel_y." + name + ".mesh" ).data() );
        system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
                  " vel_z." + name + ".mesh" ).data() );
    }
//    }



}



} // namespace LifeV


#endif //_OSEEN_H_
