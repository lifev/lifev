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

  \brief This file contains a Advection Diffusion Reaction equation solver class
*/
#ifndef _ADR_H_
#define _ADR_H_

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
#include <life/lifesolver/dataADR.hpp>
//
#include <boost/shared_ptr.hpp>

#include "life/lifefem/FESpace.hpp"

#include <list>

namespace LifeV
{
/*!
  \class ADRSolver

  This class implements a NavierStokes solver.
  The resulting linear systems are solved by GMRES on the full
  matrix ( u and p coupled ).

*/



template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class ADRSolver
//     :
//     public NavierStokesHandler<Mesh>, EpetraHandler
{

public:

    typedef DataADR<Mesh>  data_type;

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&,const Real&, const ID& );
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
    ADRSolver( const data_type&          dataType,
               FESpace<Mesh, EpetraMap>& FESpace,
               FESpace<Mesh, EpetraMap>& betaFESpace,
               BCHandler&                bcHandler,
               Epetra_Comm&              comm );

    /*!
      \param dataType
      \param localMap
      \param velocity FE space
      \param pressure FE space
    */
    ADRSolver( const data_type&      dataType,
               FESpace<Mesh, EpetraMap>& FESpace,
               FESpace<Mesh, EpetraMap>& betaFESpace,
               Epetra_Comm&              comm );

    //! virtual destructor

    virtual ~ADRSolver();

    //! Update convective term, bc treatment and solve the linearized ns system

    virtual void iterate( bchandler_raw_type& bch );

    virtual void setUp        ( const GetPot& dataFile );

    virtual void buildSystem();

    virtual void updateRHS(vector_type const& sourceVec)
    {
        M_rhsNoBC = sourceVec;
        M_rhsNoBC.GlobalAssemble();
    }

    virtual void updateSystem(Real       alpha,
                              vector_type& betaVec,
                              vector_type& sourceVec
                              );

//     void updateLinearSystem(Real       alpha,
//                             vector_type& betaVec,
//                             vector_type& sourceVec
//                             );

    void initialize( const Function& );
    void initialize( const vector_type& );

    //! returns the local solution vector

    const vector_type& solution() const {return M_sol;}

    //! returns the local residual vector

    const vector_type& residual() const {return M_residual;}

    //! reduce the local solution solution in global vectors

    void reduceSolution( Vector& u,
                         Vector& p );

    void reduceResidual( Vector& res);


    FESpace<Mesh, EpetraMap>& velFESpace()   {return M_FESpace;}


    //! Bounday Conditions
    /*const*/ bool BCset() const {return M_setBC;}
    //! set the BCs
    void setBC(BCHandler &BCh)
        {
            M_BCh = &BCh; M_setBC = true;
        }

//     BCHandler& bcHandler()
//         {
//             return *M_BCh;
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

    // compute area on a face with given flag
    Real area(const EntityFlag& flag);

    // compute flux on a face with given flag
    Real flux(const EntityFlag& flag);

    //! Postprocessing
    void postProcess(bool _writeMesh = false);

    void resetPrec() {M_resetPrec = true; M_resetStab = true;}
    // as for now resetting stabilization matrix at the same time as the preconditioner
    // void resetStab() {M_resetStab = true;}

    //Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.getRepeatedEpetra_Map(); }

    EpetraMap const& getMap() const { return M_localMap; }

    const Epetra_Comm& comm() const {return *M_comm;}

    bool isLeader() const
    {
        assert( M_comm != 0);
        return comm().MyPID() == 0;
    }

    void leaderPrint   (string const message, Real const number) const;
    void leaderPrint   (string const message) const;
    void leaderPrintMax(string const message, Real const number) const;


    void recomputeMatrix(bool const recomp){M_recomputeMatrix = recomp;}

    matrix_type& matrNoBC()
        {
            return *M_matrNoBC;
        }
    matrix_type& matrMass()
        {
            return *M_matrMass;
        }


protected:

    UInt dim() const           { return M_FESpace.dim(); }


    void solveSystem            (  matrix_ptrtype matrFull,
                                   vector_type&   rhsFull,
                                   vector_type&    sol,
                                   SolverType&     linearSolver,
                                   prec_type&      prec );

    void applyBoundaryConditions(  matrix_type&        matrix,
                                   vector_type&        rhs,
                                   bchandler_raw_type& BCh);

    void echo(std::string message);

    //! removes mean of component comp of vector x
    Real removeMean( vector_type& x );

    //private members

    //! data for NS solvers
    const data_type&               M_data;

    // FE spaces
    FESpace<Mesh, EpetraMap>&      M_FESpace;

    // beta FE spaces
    FESpace<Mesh, EpetraMap>&      M_betaFESpace;

    //! MPI communicator
    Epetra_Comm*                   M_comm;
    int                            M_me;

    //! Bondary Conditions Handler
    BCHandler*                     M_BCh;
    bool                           M_setBC;

    EpetraMap                      M_localMap;

    //! mass matrix

    matrix_ptrtype                 M_matrMass;

    //! Stiffness Matrix mu*stiff
    matrix_ptrtype                 M_matrStiff;

    //! Advection Matrix
    matrix_ptrtype                 M_matrAdv;

    //! matrix to be solved
//    matrix_ptrtype                 M_matrFull;

    //! matrix without boundary conditions

    matrix_ptrtype                 M_matrNoBC;

    //! stabilization matrix
    matrix_ptrtype                 M_matrStab;

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

     boost::shared_ptr<EpetraPreconditioner> M_prec;

    bool                           M_steady;

    //! Stabilization
    std::string                    M_stab;
    Real                           M_gammaBeta;
    bool                           M_resetStab;
    bool                           M_reuseStab;

//     details::
//    IPStabilization<Mesh, Dof>     M_ipStab;

    const Function*                M_betaFct;

    bool                           M_divBetaUv;

    //
    Real                         M_diagonalize;

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
    ElemMat                        M_elmatStiff;      // stiffness
    ElemMat                        M_elmatAdv;        // advection
    ElemMat                        M_elmatMass;       // mass
    ElemMat                        M_elmatStab;
//    ElemVec                        M_elvec;           // Elementary right hand side
    ElemVec                        M_elvec_u;         // Elementary right hand side


}; // class ADRSolver



//
// IMPLEMENTATION
//


template<typename Mesh, typename SolverType>
ADRSolver<Mesh, SolverType>::
ADRSolver( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& FESpace,
           FESpace<Mesh, EpetraMap>& betaFESpace,
           BCHandler&                BCh,
           Epetra_Comm&              comm ):
    M_data                   ( dataType ),
    M_FESpace                ( FESpace ),
    M_betaFESpace            ( betaFESpace ),
    M_comm                   ( &comm ),
    M_me                     ( M_comm->MyPID() ),
    M_BCh                    ( &BCh ),
    M_setBC                  ( true ),
    M_localMap               ( M_FESpace.map() ),
    M_matrMass               ( ),
    M_matrStiff              ( ),
//    M_matrFull               ( ),
    M_matrNoBC               ( ),
    M_rhsNoBC                ( M_localMap ),
    M_rhsFull                ( M_localMap ),
    M_sol                    ( M_localMap ),
    M_residual               ( M_localMap ),
    M_linearSolver           ( ),
    M_prec                   ( ),
    M_resetStab              ( true ),
    M_reuseStab              ( true ),
//     M_ipStab                 ( M_FESpace.mesh(),
//                                M_FESpace.dof(), M_FESpace.refFE(),
//                                M_FESpace.feBd(), M_FESpace.qr(),
//                                0., 0., 0.,
//                                M_data.diffusivity() ),
    M_betaFct                ( 0 ),
    M_count                  ( 0 ),
    M_verbose                ( M_me == 0),
    M_updated                ( false ),
    M_reusePrec              ( true ),
    M_maxIterForReuse        ( -1 ),
    M_resetPrec              ( true ),
    M_maxIterSolver          ( -1 ),
    M_recomputeMatrix        ( false ),
    M_elmatStiff             ( M_FESpace.fe().nbNode, 1, 1 ),
    M_elmatMass              ( M_FESpace.fe().nbNode, 1, 1 ),
    M_elmatAdv               ( M_FESpace.fe().nbNode, 1, 1 ),
    M_elmatStab              ( M_FESpace.fe().nbNode, 1, 1 ),
//    M_elvec                  ( M_FESpace.fe().nbNode, nDimensions ),
    M_elvec_u                ( M_betaFESpace.fe().nbNode, nDimensions ) //SQ: from M_FESpace to M_betaFESpace
{
}



template<typename Mesh, typename SolverType>
ADRSolver<Mesh, SolverType>::
ADRSolver( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& FESpace,
           FESpace<Mesh, EpetraMap>& betaFESpace,
           Epetra_Comm&              comm ):
    M_data                   ( dataType ),
    M_FESpace                ( FESpace ),
    M_betaFESpace            ( betaFESpace ),
    M_comm                   ( &comm ),
    M_me                     ( M_comm->MyPID() ),
    M_setBC                  ( false ),
    M_localMap               ( M_FESpace.map() ),
    M_matrMass               ( ),
    M_matrStiff              ( ),
    M_matrNoBC               ( ),
    M_rhsNoBC                ( M_localMap ),
    M_rhsFull                ( M_localMap ),
    M_sol                    ( M_localMap ),
    M_residual               ( M_localMap ),
    M_linearSolver           ( ),
    M_prec                   ( ),
    M_resetStab              ( true ),
    M_reuseStab              ( true ),
    M_betaFct                ( 0 ),
    M_count                  ( 0 ),
    M_verbose                ( M_me == 0),
    M_updated                ( false ),
    M_reusePrec              ( true ),
    M_maxIterForReuse        ( -1 ),
    M_resetPrec              ( true ),
    M_maxIterSolver          ( -1 ),
    M_recomputeMatrix        ( false ),
    M_elmatStiff             ( M_FESpace.fe().nbNode, 1, 1 ),
    M_elmatMass              ( M_FESpace.fe().nbNode, 1, 1 ),
    M_elmatAdv               ( M_FESpace.fe().nbNode, 1, 1 ),
    M_elmatStab              ( M_FESpace.fe().nbNode, 1, 1 ),
    M_elvec_u                ( M_betaFESpace.fe().nbNode, nDimensions ) //SQ: from M_FESpace to M_betaFESpace
{
}



template<typename Mesh, typename SolverType>
ADRSolver<Mesh, SolverType>::
~ADRSolver()
{

}

template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
leaderPrint(string const message, Real const number) const
{
  if ( isLeader() )
    std::cout << message << number << std::endl;

}

template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
leaderPrint(string const message) const
{
  if ( isLeader() )
    std::cout << message << std::flush;

}

template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
leaderPrintMax(string const message, Real const number) const
{
  Real num(number);
  Real globalMax;
  M_comm->MaxAll(&num, &globalMax, 1);

  leaderPrint( message , globalMax );

}



template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::setUp( const GetPot& dataFile )
{
    M_steady      = dataFile( "adr/miscellaneous/steady",        0  );
    M_stab        = dataFile( "adr/stab/type",                   "ip");
    M_gammaBeta   = dataFile( "adr/stab/gammaBeta",              0. );
    M_reuseStab   = dataFile( "adr/stab/reuse",                  true);
    M_diagonalize = dataFile( "adr/space_discretization/diagonalize",  0. );


    M_linearSolver.setDataFromGetPot( dataFile, "adr/solver" );

    M_maxIterSolver   = dataFile( "adr/solver/max_iter", -1);
    M_reusePrec       = dataFile( "adr/prec/reuse", true);
    M_maxIterForReuse = dataFile( "adr/prec/max_iter_reuse", M_maxIterSolver*8/10);

    std::string precType = dataFile( "adr/prec/prectype", "Ifpack");

    M_prec.reset( PRECFactory::instance().createObject( precType ) );
    ASSERT(M_prec.get() != 0, "AdvectionDiffusionSolver : Preconditioner not set");

    M_prec->setDataFromGetPot( dataFile, "adr/prec" );
}

template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::buildSystem()
{


    M_matrMass.reset  ( new matrix_type(M_localMap) );
    M_matrStiff.reset( new matrix_type(M_localMap) );

//    M_comm->Barrier();

    leaderPrint("  adr-  Computing constant matrices ...          ");

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
    UInt nbCompU = 1;

    // Elementary computation and matrix assembling
    // Loop on elements

    UInt velTotalDof   = M_FESpace.dof().numTotalDof();


    chrono.start();

    for ( UInt iVol = 1; iVol <= M_FESpace.mesh()->numElements(); iVol++ )
    {
        chronoDer.start();
        M_FESpace.fe().updateFirstDeriv( M_FESpace.mesh()->element( iVol ) );

        chronoDer.stop();

        chronoZero.start();
        M_elmatStiff.zero();
        M_elmatMass.zero();
        chronoZero.stop();


        // stiffness strain
        chronoStiff.start();
        //stiff_strain( 2.0*M_data.viscosity(), M_elmatStiff, M_FESpace.fe() );
        stiff( M_data.diffusivity(), M_elmatStiff,  M_FESpace.fe(), 0, 0 );
        //stiff_div( 0.5*M_FESpace.fe().diameter(), M_elmatStiff, M_FESpace.fe() );
        chronoStiff.stop();

        // mass
        if ( !M_steady )
        {
            chronoMass.start();
            mass( 1., M_elmatMass, M_FESpace.fe(), 0, 0);
            chronoMass.stop();
        }

        // stiffness

        chronoStiffAssemble.start();
        assembleMatrix( *M_matrStiff,
                        M_elmatStiff,
                        M_FESpace.fe(),
                        M_FESpace.fe(),
                        M_FESpace.dof(),
                        M_FESpace.dof(),
                        0, 0,
                        0, 0);
        chronoStiffAssemble.stop();

        if ( !M_steady )
        {
            chronoMassAssemble.start();
            assembleMatrix( *M_matrMass,
                            M_elmatMass,
                            M_FESpace.fe(),
                            M_FESpace.fe(),
                            M_FESpace.dof(),
                            M_FESpace.dof(),
                            0, 0, 0, 0);
            chronoMassAssemble.stop();
        }

    }



    M_comm->Barrier();

    chrono.stop();
    leaderPrintMax( "done in " , chrono.diff());


    leaderPrint( "  adr-  Finalizing the matrices ...              ");
    chrono.start();

    M_matrStiff->GlobalAssemble();
    M_matrMass->GlobalAssemble();

//     M_matrStiff->spy("stiff");
//     M_matrMass->spy("mass");

    chrono.stop();
    leaderPrintMax("done in " , chrono.diff() );;

    if (false)
        std::cout << "partial times:  \n"
                  << " Der            " << chronoDer.diff_cumul() << " s.\n"
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
void ADRSolver<Mesh, SolverType>::
initialize( const Function& u0 )
{
     vector_type u(M_FESpace.map());
     M_FESpace.interpolate(u0, u, 0.0);
     M_sol = u;

}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
initialize( const vector_type& u0 )
{
    M_sol = u0;
}

template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
updateSystem( Real       alpha,
              vector_type& betaVec,
              vector_type& sourceVec
              )
{

    Chrono chrono;

    leaderPrint("  adr-  Updating mass term on right hand side... ");

    chrono.start();

    UInt velTotalDof   = M_FESpace.dof().numTotalDof();

    // Right hand side for the velocity at time

    updateRHS(sourceVec);

    chrono.stop();

    leaderPrintMax("done in ", chrono.diff());


    M_updated = false;

//

    if (M_recomputeMatrix)
        buildSystem();

    leaderPrint( "  adr-  Copying the matrices ...                 ");

    chrono.start();

     if (M_matrNoBC)
         M_matrNoBC.reset(new matrix_type(M_localMap, M_matrNoBC->getMeanNumEntries() ));
     else
         M_matrNoBC.reset(new matrix_type(M_localMap, M_matrStiff->getMeanNumEntries() ));

     M_matrStab.reset( new matrix_type(M_localMap) );

    if ( M_data.diffusivity() != 0. )
        *M_matrNoBC += *M_matrStiff;

    if (alpha != 0. )
    {
        *M_matrNoBC += *M_matrMass*alpha;
    }


    chrono.stop();
    leaderPrintMax( "done in " , chrono.diff() );

//    UInt nbCompU       = nDimensions;

    //! managing the convective term

    Real normInf;
    betaVec.NormInf(&normInf);

    if (normInf != 0.)
    {

        leaderPrint("  adr-  Updating the convective terms ...        ");

        // vector with repeated nodes over the processors

        vector_type betaVecRep(betaVec, Repeated );
        int uDim = betaVec.size()/nDimensions;
        chrono.start();

        for ( UInt iVol = 1; iVol <= M_FESpace.mesh()->numElements(); ++iVol )
        {

            M_FESpace.fe().updateFirstDeriv( M_FESpace.mesh()->element( iVol ) ); //as updateFirstDer
            M_betaFESpace.fe().updateFirstDeriv( M_FESpace.mesh()->element( iVol ) ); //as updateFirstDer

            M_elmatAdv.zero();

            UInt eleID = M_betaFESpace.fe().currentLocalId();
            // Non linear term, Semi-implicit approach
            // M_elvec contains the velocity values in the nodes
            for ( UInt iNode = 0 ; iNode < ( UInt ) M_betaFESpace.fe().nbNode ; iNode++ )
            {
	         UInt  iloc = M_betaFESpace.fe().patternFirst( iNode );
                for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
                {
		    UInt ig = M_betaFESpace.dof().localToGlobal( eleID, iloc + 1 ) + iComp*uDim;
                    M_elvec_u.vec()[ iloc + iComp*M_betaFESpace.fe().nbNode ] = betaVecRep[ig]; // BASEINDEX + 1
                }
            }


            // compute local convective term and assembling
            for(UInt iComp=0; iComp<nDimensions; iComp++)
            	grad( iComp, M_elvec_u, M_elmatAdv, M_FESpace.fe(), M_FESpace.fe(), M_betaFESpace.fe());




            assembleMatrix( *M_matrNoBC,
                            M_elmatAdv,
                            M_FESpace.fe(),
                            M_FESpace.fe(),
                            M_FESpace.dof(),
                            M_FESpace.dof(),
                            0, 0,
                            0, 0
                            );


            // Streamline Diffusion ( only stab for now )


            if (M_stab == "sd")
            {
                M_elmatStab.zero();
                Real VLoc_infty = 0.;
                Real VLoc_mean  = 0.;
                Real VLoc_c     = 0.;

                for ( UInt ih_c = 0 ; ih_c < ( UInt ) this->M_betaFESpace.fe().nbNode ; ih_c++ )
                {
                    UInt iloc = this->M_betaFESpace.fe().patternFirst( ih_c );
                    for ( UInt iComp = 0; iComp < nDimensions; ++iComp)
                    {
                        UInt ig = M_betaFESpace.dof().localToGlobal( eleID, iloc + 1 ) + iComp*uDim;
                        M_elvec_u.vec()[ iloc + iComp * this->M_betaFESpace.fe().nbNode ] = betaVecRep[ ig ];
                        VLoc_c += betaVecRep[ ig ] * betaVecRep[ ig ];
                    }

                    VLoc_c     = sqrt( VLoc_c );
                    VLoc_mean += VLoc_c;

                    if ( VLoc_c > VLoc_infty )
                        VLoc_infty = VLoc_c;
                }

                VLoc_mean = VLoc_mean / this->M_betaFESpace.fe().nbNode;

                Real coef_stab, Pe_loc = 0;
                coef_stab=M_gammaBeta*this->M_betaFESpace.fe().diameter()*VLoc_infty; // Alessandro - method

                stiff_sd( coef_stab / ( VLoc_mean*VLoc_mean ), M_elvec_u, M_elmatStab, this->M_FESpace.fe(), this->M_betaFESpace.fe() );

                assembleMatrix( *M_matrStab,
                                M_elmatStab,
                                M_FESpace.fe(),
                                M_FESpace.fe(),
                                M_FESpace.dof(),
                                M_FESpace.dof(),
                                0, 0, 0, 0 );
				}
        }

//TODO: check both the 2D and 3D implementation of ip-stabilization (on the Laplacian test doesn't work properly)


        if (M_stab == "ip")
        {
//	      leaderPrint("   adr- IP stab");
            if ( M_resetStab )
            {
                const UInt nDof = M_betaFESpace.dof().numTotalDof();

                CurrentFE fe1(M_FESpace.refFE(),
                              getGeoMap(*M_FESpace.mesh()),
                              M_FESpace.qr());
                CurrentFE fe2(M_FESpace.refFE(),
                              getGeoMap(*M_FESpace.mesh()),
                              M_FESpace.qr());
                CurrentFE fe3(M_betaFESpace.refFE(),
                              getGeoMap(*M_betaFESpace.mesh()),
                              M_betaFESpace.qr());


#ifdef TWODIM
                typedef ID ( *ETOP )( ID const localFace, ID const point );
                ETOP  eToP;
                switch( M_FESpace.fe().refFE.type )
                {
                    case FE_P1_2D:
                    	eToP = LinearTriangle::eToP;
                        break;
                    case FE_P2_2D:
                    	eToP = QuadraticTriangle::eToP;
                        break;
                    case FE_Q1_2D:
                        eToP = LinearQuad::eToP;
                        break;
                    case FE_Q2_2D:
                        eToP = QuadraticQuad::eToP;
                        break;
                    default:
                    	eToP=0;
                        ERROR_MSG( "This refFE is not allowed with IP stabilization" );
                        break;
                }
                for ( UInt iEdge = M_FESpace.mesh()->numBEdges() + 1; iEdge <= M_FESpace.mesh()->numEdges();
                      ++iEdge )
                {
                    const UInt iElAd1 = M_FESpace.mesh()->edge( iEdge ).ad_first();
                    const UInt iElAd2 = M_FESpace.mesh()->edge( iEdge ).ad_second();

                    if ( iElAd1 == iElAd2 || iElAd1 == 0 || iElAd2 == 0)
                    {
                        continue;
                    }

                    M_elmatStab.zero();

                    M_betaFESpace.feBd().updateMeas( M_betaFESpace.mesh()->edge( iEdge ) );
                    const Real hK2  = std::pow(M_betaFESpace.feBd().measure(), 2.);

                    M_betaFESpace.feBd().updateMeasNormal( M_betaFESpace.mesh()->edge( iEdge ) );
                    KNM<Real>& normal = M_betaFESpace.feBd().normal;

                    fe1.updateFirstDeriv( M_FESpace.mesh()->element( iElAd1 ) );
                    fe2.updateFirstDeriv( M_FESpace.mesh()->element( iElAd2 ) );

                   ElemVec beta(M_betaFESpace.feBd().nbNode, nDimensions);

                    // first, get the local trace of the velocity into beta
                    // local id of the face in its adjacent element

                    UInt iEdEl = M_betaFESpace.mesh()->edge( iEdge ).pos_first();
                    for ( int iNode = 0; iNode < M_betaFESpace.feBd().nbNode; ++iNode )
                    {
                        UInt iloc = eToP( iEdEl, iNode+1 );
                        for ( int iCoor = 0; iCoor < fe1.nbCoor; ++iCoor )
                        {
                            UInt ig = M_betaFESpace.dof().localToGlobal( iElAd1, iloc + 1 ) - 1 +iCoor*nDof;
                            if (betaVecRep.BlockMap().LID(ig + 1) >= 0)
                                beta.vec()[ iCoor*M_betaFESpace.feBd().nbNode + iNode ] = betaVecRep( ig + 1); // BASEINDEX + 1
                        }
                    }

#elif defined THREEDIM
                typedef ID ( *FTOP )( ID const localFace, ID const point );
                 FTOP  fToP;
                switch( M_FESpace.fe().refFE.type )
                {
                case FE_P1_3D:
                case FE_P1bubble_3D:
                	fToP = LinearTetra::fToP;
                	break;
                case FE_P2_3D:
                	fToP = QuadraticTetra::fToP;
                	break;
                case FE_Q1_3D:
					fToP = LinearHexa::fToP;
					break;
                case FE_Q2_3D:
                	fToP = QuadraticHexa::fToP;
                	break;
				default:
					fToP = 0;
					ERROR_MSG( "This refFE is not allowed with IP stabilisation" );
					break;
                }


                for ( UInt iFace = M_FESpace.mesh()->numBFaces() + 1; iFace <= M_FESpace.mesh()->numFaces();
                      ++iFace )
                {

                    const UInt iElAd1 = M_FESpace.mesh()->face( iFace ).ad_first();
                    const UInt iElAd2 = M_FESpace.mesh()->face( iFace ).ad_second();

                    if ( iElAd1 == iElAd2 || iElAd1 == 0 || iElAd2 == 0)
                    {
                        continue;
                    }

                    M_elmatStab.zero();

                    M_betaFESpace.feBd().updateMeas( M_betaFESpace.mesh()->face( iFace ) );
                    const Real hK2  = std::pow(M_betaFESpace.feBd().measure(), 2.);

                    M_betaFESpace.feBd().updateMeasNormal( M_betaFESpace.mesh()->face( iFace ) );
                    KNM<Real>& normal = M_betaFESpace.feBd().normal;

                    fe1.updateFirstDeriv( M_FESpace.mesh()->element( iElAd1 ) );
                    fe2.updateFirstDeriv( M_FESpace.mesh()->element( iElAd2 ) );

                    Real bn   = 0;
                    Real bmax = 0;

		    // Old version, removed by SQ
		    /*
                    ElemVec beta( M_betaFESpace.feBd().nbNode, nDimensions);

                    // first, get the local trace of the velocity into beta
                    // local id of the face in its adjacent element

                    UInt iFaEl = M_betaFESpace.mesh()->face( iFace ).pos_first();
                    for ( int iNode = 0; iNode < M_betaFESpace.feBd().nbNode; ++iNode )
                    {
                        UInt iloc = fToP( iFaEl, iNode+1 );
                        for ( int iCoor = 0; iCoor < fe1.nbCoor; ++iCoor )
                        {
                            UInt ig = M_betaFESpace.dof().localToGlobal( iElAd1, iloc + 1 ) - 1 +iCoor*nDof;
                            if (betaVecRep.BlockMap().LID(ig + 1) >= 0)
			      {
				Real value( betaVecRep( ig + 1) );
                                beta.vec()[ iCoor*M_betaFESpace.feBd().nbNode + iNode ] = value; // BASEINDEX + 1
			      };
                        }
		    }*/
		    
		    // New version, added by SQ : beta on the domain, not the boundary!
		    // See elemOper.cpp for justification of this usage
                    ElemVec beta( M_betaFESpace.fe().nbNode, nDimensions);

                    for ( int iNode = 0; iNode < M_betaFESpace.fe().nbNode; ++iNode )
                    {
		        UInt  iloc = M_betaFESpace.fe().patternFirst( iNode );
                        for ( int iCoor = 0; iCoor < fe1.nbCoor; ++iCoor )
                        {
                            UInt ig = M_betaFESpace.dof().localToGlobal( iElAd1, iloc + 1 ) + iCoor*nDof;
			    beta.vec()[ iloc + iCoor*M_betaFESpace.fe().nbNode ] = betaVecRep[ig]; // BASEINDEX + 1

                        }
		    }

                    // second, calculate its max norm
                    for ( int l = 0; l < int( M_betaFESpace.fe().nbCoor*M_betaFESpace.fe().nbNode ); ++l ) // SQ: feBd->fe
                    {
		      if ( bmax < fabs( beta.vec()[ l ] ) )
			bmax = fabs( beta.vec()[ l ] );
                    }

                     UInt iFaEl = M_betaFESpace.mesh()->face( iFace ).pos_first();
                     for ( int iNode = 0; iNode < M_betaFESpace.feBd().nbNode; ++iNode )
                     {
                         UInt iloc = fToP( iFaEl, iNode + 1 );
                         for ( int iCoor = 0; iCoor < nDimensions; ++iCoor )
                         {
                             UInt ig = M_betaFESpace.dof().localToGlobal( iElAd1, iloc + 1 ) - 1 + iCoor*nDof;
                             if (betaVecRep.BlockMap().LID(ig + 1) >= 0)
                                 bn += normal(iNode, iCoor)*betaVecRep( ig + 1 );
                         }
                     }

                    Real coeffBeta = hK2*M_gammaBeta*abs(bn);


#endif
		    //                    Real coeffBeta = M_gammaBeta;
		    //                     std::cout << coeffBeta << std::endl;

                    ipstab_bagrad( coeffBeta,
                                   M_elmatStab,
                                   fe1,
                                   fe2,
                                   fe3,
                                   beta,
                                   M_FESpace.feBd(),
                                   0, 0);

                    assembleMatrix( *M_matrStab,
                                    M_elmatStab,
                                    fe1,
                                    fe2,
                                    M_FESpace.dof(),
                                    M_FESpace.dof(),
                                    0, 0, 0, 0 );

                    M_resetStab = true;
                }
            }
        }

        chrono.stop();
        leaderPrintMax( "done in " , chrono.diff() );
    }


    M_updated = true;

}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::iterate( bchandler_raw_type& bch )
{

    Chrono chrono;


    // matrix and vector assembling communication

    leaderPrint("  adr-  Finalizing the matrix and vectors ...    ");

    chrono.start();

    M_matrNoBC->GlobalAssemble();
    M_matrStab->GlobalAssemble();

    chrono.stop();
    leaderPrintMax("done in ", chrono.diff() );

    //

    leaderPrint("  adr-  setting up the full matrix        ...    ");

    chrono.start();

    matrix_ptrtype matrFull( new matrix_type( M_localMap, M_matrNoBC->getMeanNumEntries()));

    *matrFull += *M_matrNoBC;
    *matrFull += *M_matrStab;

    chrono.stop();
    leaderPrintMax("done in ", chrono.diff() );

    //
    leaderPrint("  adr-  setting up the full rhs           ...    ");

    chrono.start();


//    M_matrStab->spy("stab");
    vector_type    rhsFull = M_rhsNoBC;

    chrono.stop();
    leaderPrintMax("done in ", chrono.diff() );

    // boundary conditions update

    leaderPrint("  adr-  Applying boundary conditions ...         ");
    chrono.start();

    applyBoundaryConditions( *matrFull, rhsFull, bch);

    M_comm->Barrier();
    leaderPrintMax("done in " , chrono.diff());

    //

    leaderPrint("  adr-  Finalizing the full matrix    ...        ");
    chrono.start();


    matrFull->GlobalAssemble();

    chrono.stop();

    M_comm->Barrier();

    leaderPrintMax("done in " , chrono.diff());
    // solving the system

    solveSystem( matrFull, rhsFull, M_sol, M_linearSolver, M_prec);

    M_residual  = M_rhsNoBC;
    M_residual -= *M_matrNoBC*M_sol;

} // iterate()



template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::solveSystem( matrix_ptrtype  matrFull,
                                           vector_type&    rhsFull,
                                           vector_type&    sol,
                                           SolverType&     linearSolver,
                                           prec_type&      prec)
{
    Chrono chrono;

    leaderPrint("  adr-  Setting up the solver ...                ");

    chrono.start();
    linearSolver.setMatrix(*matrFull);
    chrono.stop();

    leaderPrintMax("done in " , chrono.diff());

    // overlapping schwarz preconditioner

    if ( !M_reusePrec || M_resetPrec || !prec->set() )
    {
        chrono.start();

        leaderPrint("  adr-  Computing the precond ...                ");

        prec->buildPreconditioner(matrFull);

        Real condest = prec->Condest();

        linearSolver.setPreconditioner(prec);

        chrono.stop();
        leaderPrintMax( "done in " , chrono.diff() );
	leaderPrint("  adr-       Estimated condition number = " , condest );

    }
    else
    {
        leaderPrint("  adr-  Reusing  precond ...                \n");
    }


    chrono.start();

    leaderPrint("  adr-  Solving system ...                       ");

    int numIter = linearSolver.solve(sol, rhsFull);

    if (numIter > M_maxIterSolver)
    {
        chrono.start();

	leaderPrint("  adr- Iterative solver failed, numiter = " , numIter);
        leaderPrint("     maxIterSolver = " , M_maxIterSolver );
	leaderPrint("     recomputing the precond ...            ");

        prec->buildPreconditioner(matrFull);

        Real condest = prec->Condest();

        linearSolver.setPreconditioner(prec);

        chrono.stop();
        leaderPrintMax( "done in " , chrono.diff() );
        leaderPrint("  adr-       Estimated condition number = " , condest );

        numIter = linearSolver.solve(sol, rhsFull);

        if (numIter > M_maxIterSolver && M_verbose)
            std::cout << "  adr- ERROR: Iterative solver failed again.\n" <<  std::flush;

    }

    M_resetPrec = (numIter > M_maxIterForReuse);

    leaderPrintMax( "done in " , chrono.diff() );

    std::string statusMessage;

    statusMessage = linearSolver.printStatus();

    leaderPrint("  adr-  " + statusMessage);

    M_comm->Barrier();

}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::reduceSolution( Vector& u,
                                              Vector& p )
{
    vector_type vel(M_sol, 0);

    if (M_verbose)
    {
        for ( UInt iDof = 0; iDof < nDimensions*dim(); ++iDof )
        {
            u[ iDof ] = vel[ iDof + 1 ]; // BASEINDEX + 1
        }
    }

}

template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::reduceResidual( Vector& res )
{
    vector_type vel(M_residual, 0);

    if (M_verbose)
    {
        for ( UInt iDof = 0; iDof < nDimensions*dim(); ++iDof )
        {
            res[ iDof ] = vel[ iDof + 1 ]; // BASEINDEX + 1
        }

    }
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::applyBoundaryConditions( matrix_type&        matrix,
                                                       vector_type&        rhs,
                                                       bchandler_raw_type& BCh )
{

    // M_rhsFull = M_rhsNoBC;

    // BC manage for the velocity
    Chrono chrono;

    if ( !BCh.bdUpdateDone() )
    {

        leaderPrint( "\n     - Updating the BC ... ");
        chrono.start();
        BCh.bdUpdate( *M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );
        chrono.stop();
        leaderPrintMax( "done in " , chrono.diff() );
        leaderPrint( "\n");
    }

    //    vector_type rhsFull(rhs, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately
    vector_type rhsFull(rhs, Unique); // ignoring non-local entries, Otherwise they are summed up lately


    leaderPrint( "\n     - Managing the BC ... ");
    chrono.start();
    bcManage( matrix, rhsFull, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 1.,
              M_data.getTime() );
    chrono.stop();
    leaderPrintMax( "done in " , chrono.diff() );


    rhs = rhsFull;


} // applyBoundaryCondition


//! Computes the flux on a given part of the boundary
template<typename Mesh, typename SolverType> Real
ADRSolver<Mesh, SolverType>::flux(const EntityFlag& flag) {

  Real flux = 0.0;

  PhysVectUnknown<Vector> u(nDimensions*dim());

  reduceSolution(u);

  if (M_verbose)
  {
	  typedef  typename Mesh::ElementShape GeoShape;
	  typedef typename GeoShape::GeoBShape GeoBShape;

	  // Some useful local variables, to save some typing
	  UInt nDofPerVert = M_FESpace.refFE().nbDofPerVertex; // number of Dof per vertices
	  UInt nDofPerEdge = M_FESpace.refFE().nbDofPerEdge;   // number of Dof per edges
	  UInt nDofPerFace = M_FESpace.refFE().nbDofPerFace;   // number of Dof per faces

	  UInt nBElementV = GeoBShape::numVertices; // Number of face's vertices
	  UInt nBElementE = GeoBShape::numEdges;    // Number of face's edges

	  UInt nElemV = GeoShape::numVertices; // Number of element's vertices
	  UInt nElemE = GeoShape::numEdges;    // Number of element's edges

	  UInt nDofBElV = nDofPerVert * nBElementV; // number of vertex's Dof on a face
	  UInt nDofBElE = nDofPerEdge * nBElementE; // number of edge's Dof on a face

#ifdef TWODIM
    UInt nDofBEl = nDofBElV + nDofBElE; // number of total Dof on a boundary element
#elif defined THREEDIM
    UInt nDofBEl = nDofBElV + nDofBElE + nDofPerFace; // number of total Dof on a boundary element
#endif

	  UInt nDofElemV = nElemV*nDofPerVert; // number of vertex's Dof on a Element
	  UInt nDofElemE = nElemE*nDofPerEdge; // number of edge's Dof on a Element

	  UInt bdnF  = M_data.mesh()->numBElements();    // number of faces on boundary

	  std::list<std::pair<ID, SimpleVect<ID> > > faces;
	  ID ibF;
	  UInt iElAd, iVeEl, iBElEl, iEdEl;
	  ID lDof, gDof, numTotalDof=M_FESpace.dof().numTotalDof();

	  EntityFlag marker;
	  typedef std::list<pair<ID, SimpleVect<ID> > >::iterator Iterator;

	  //
	  // Loop on boundary faces: List of boundary faces
	  // with marker = flag
	  //
	  for (ID i=1 ; i<=bdnF; ++i) {
	    marker = M_data.mesh()->bElement(i).marker();
	    if ( marker == flag  ) {
	      faces.push_front(make_pair(i,SimpleVect<ID>(nDofBEl)));
	    }
	  }

	  //
	  // Loop on faces: building the local to global vector
	  // for these boundary faces
	  //
	  for (Iterator j=faces.begin(); j != faces.end(); ++j) {

	    ibF = j->first;

	    iElAd = M_data.mesh()->bElement(ibF).ad_first();  // id of the element adjacent to the face
	    iBElEl = M_data.mesh()->bElement(ibF).pos_first(); // local id of the boundary Element in its adjacent element

	    // Vertex based Dof
	    if ( nDofPerVert ) {

	      // loop on face vertices
	      for (ID iVeBEl=1; iVeBEl<=nBElementV; ++iVeBEl){


#ifdef TWODIM
				iVeEl = GeoShape::eToP( iBElEl, iVeBEl ); // local vertex number (in element)
#elif defined THREEDIM
                iVeEl = GeoShape::fToP( iBElEl, iVeBEl ); // local vertex number (in element)
#endif

	      	// Loop number of Dof per vertex
	      	for (ID l=1; l<=nDofPerVert; ++l) {
	      		lDof =   (iVeBEl-1) * nDofPerVert + l ; // local Dof j-esimo grado di liberta' su una faccia
	      		gDof =  M_FESpace.dof().localToGlobal( iElAd, (iVeEl-1)*nDofPerVert + l); // global Dof
	      		j->second( lDof ) =  gDof; // local to global on this face
	      	}
	      }
	    }

	    // Edge based Dof
	    if (nDofPerEdge) {

	      // loop on boundary element edges
	      for (ID iEdBEl=1; iEdBEl<=nBElementE; ++iEdBEl) {
#ifdef TWODIM
              iEdEl = iBElEl; // local edge number (in element)
#elif defined THREEDIM
                iEdEl = GeoShape::fToE( iBElEl, iEdBEl ).first; // local edge number (in element)
#endif
	      		// Loop number of Dof per edge
	      	for (ID l=1; l<=nDofPerEdge; ++l) {

	      		lDof =  nDofBElV + (iEdBEl-1) * nDofPerEdge + l ; // local Dof sono messi dopo gli lDof dei vertici
	      		gDof =  M_FESpace.dof().localToGlobal( iElAd, nDofElemV + (iEdEl-1)*nDofPerEdge + l); // global Dof
	      		j->second( lDof ) =  gDof; // local to global on this face
	      	}
	      }
	    }
#ifndef TWODIM
	    // Face based Dof
	    if (nDofPerFace) {

	      // Loop on number of Dof per face
	      for (ID l=1; l<=nDofPerFace; ++l) {
	      	lDof = nDofBElE + nDofBElV + l; // local Dof sono messi dopo gli lDof dei vertici e dopo quelli degli spigoli
	      	gDof = M_FESpace.dof().localToGlobal( iElAd, nDofElemE + nDofElemV + (iBElEl-1)*nDofPerFace + l); // global Dof
	      	j->second( lDof ) =  gDof; // local to global on this face
	      }
	    }
#endif
	  }

	  // Number of velocity components
	  UInt nc_u=nDimensions;

	  // Nodal values of the velocity in the current face
	  std::vector<Real> u_local(nc_u*nDofBEl);

	  // Loop on faces
	  for (Iterator j=faces.begin(); j != faces.end(); ++j) {

	    // Extracting nodal values of the velocity in the current face
	    for (UInt ic =0; ic<nc_u; ++ic) {
	      for (ID l=1; l<=nDofBEl; ++l) {
	      	gDof = j->second(l);
	      	u_local[ic*nDofBEl+l-1] = u(ic*numTotalDof+gDof-1);
	      }
	    }

	    // Updating quadrature data on the current face
	    M_FESpace.feBd().updateMeasNormalQuadPt(M_data.mesh()->bElement(j->first));

	    // Quadrature formula
	    // Loop on quadrature points
	    for(int iq=0; iq< M_FESpace.feBd().nbQuadPt; ++iq) {

	      // Dot product
	      // Loop on components
	      for (UInt ic =0; ic<nc_u; ++ic) {

	      	// Interpolation
	      	// Loop on local dof
	      	for (ID l=1; l<=nDofBEl; ++l)
	      		flux += M_FESpace.feBd().weightMeas(iq)
	      							* u_local[ic*nDofBEl+l-1]
	      		          * M_FESpace.feBd().phi(int(l-1),iq)
	      		          * M_FESpace.feBd().normal(int(ic),iq);
	      	}
	    	}
	  	}
  }

  return flux;
}



//! Computes the area on a given part of the boundary
template<typename Mesh, typename SolverType> Real
ADRSolver<Mesh, SolverType>::area(const EntityFlag& flag) {

  EpetraVector area( EpetraMap( M_comm->NumProc(), 1, &M_me, 0, *M_comm) );
  area *= 0.0;

	UInt bdnF  = M_data.mesh()->numBElements();    // number of faces on boundary

  std::list<ID> faces;
  typedef std::list<ID>::iterator Iterator;

  EntityFlag marker;
  //
  // Loop on boundary faces: List of boundary faces
  // with marker = flag
  //
  for (ID i=1 ; i<=bdnF; ++i) {
    marker = M_data.mesh()->bElement(i).marker();
    if ( marker == flag  ) {
      faces.push_front(i);
    }
  }

  //
  // Loop on processor faces
  //
  for (Iterator j=faces.begin(); j != faces.end(); ++j) {

    M_FESpace.feBd().updateMeas( M_data.mesh()->bElement( *j ) );  // updating finite element information

    area[M_me] += M_FESpace.feBd().measure();

  }

  Real total_area(0.);
  EpetraVector local_area( area, M_me );
  for( int i=0; i<M_comm->NumProc(); ++i )
  	total_area += local_area[i];

//  area.MeanValue( &total_area );

	return total_area;
}



// Postprocessing
template <typename Mesh, typename SolverType>
void
ADRSolver<Mesh, SolverType>::postProcess(bool /*_writeMesh*/)
{
    std::ostringstream index;
    std::ostringstream indexMe;
    std::string name;
    std::string me;


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



    vector_type concentration(M_sol, Repeated);


//         if ( fmod( float( M_count ), float( M_data.verbose() ) ) == 0.0 )
//         {
    if (M_me == 0)
        std::cout << "  ADR-  Post-processing " << std::flush;

    index << std::setfill('0') << std::setw(3);
    index << ( M_count / M_data.verbose() );
    name = index.str();

//     PhysVectUnknown<Vector> cc(dim());

//     reduceSolution(cc);

//     if (M_me == 0)
//     {
        // postprocess data file for medit
//         wr_medit_ascii_scalar( "vel_x." + name + ".bb", u.giveVec(),
//                                M_data.mesh()->numGlobalVertices() );
//         wr_medit_ascii_scalar( "vel_y." + name + ".bb", u.giveVec() + this->dim(),
//                                M_data.mesh()->numGlobalVertices() );
//         wr_medit_ascii_scalar( "vel_z." + name + ".bb", u.giveVec()+2*this->dim(),
//                                M_data.mesh()->numGlobalVertices() );
//         wr_medit_ascii_scalar( "press." + name + ".bb", p.giveVec(),
//                                p.size() );

    	// Real dt = M_data.getTimeStep();


    writeMesh("cc." + me + "." + name + ".mesh", *M_FESpace.mesh());

    meditSolutionWriter( "cc." + me + "." + name + ".bb",
                         *M_FESpace.mesh(), concentration, 0);

//     }
//    }

    M_count++;


}



} // namespace LifeV


#endif //_OSEEN_H_
