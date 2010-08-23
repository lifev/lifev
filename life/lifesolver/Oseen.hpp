/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Gilles Fourestey <gilles.fourestey@imag.fr>
            Simone Deparis   <simone.deparis@epfl.ch>
      Date: 2007-06-01

 Copyright (C) 2008 EPFL

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

#ifndef _OSEEN_H_
#define _OSEEN_H_

#include <life/lifefilters/medit_wrtrs.hpp>

#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/EpetraPreconditioner.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/EpetraMap.hpp>

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>
//
#include <life/lifecore/chrono.hpp>

#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/sobolevNorms.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/postProc.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifesolver/nsipterms.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
//#include <life/lifesolver/dataMesh.hpp>

//
#include <boost/shared_ptr.hpp>


#include <list>

namespace LifeV
{
/*!
  This file contains a Oseen equation solver class.
  The resulting linear systems are solved by GMRES on the full
  matrix ( u and p coupled ).

*/



template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class Oseen
{

public:

    typedef Mesh                                  mesh_type;
    typedef DataNavierStokes                      data_type;


//     typedef Real ( *Function ) ( const Real&, const Real&, const Real&,
//                                  const Real&, const ID& );

    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> Function;


    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> source_type;


    typedef BCHandler                             bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type> bchandler_type;

//    typedef SolverType                    solver_type;
    typedef typename SolverType::matrix_type      matrix_type;
    typedef boost::shared_ptr<matrix_type>        matrix_ptrtype;
    typedef typename SolverType::vector_type      vector_type;
    typedef boost::shared_ptr<vector_type>        vector_ptrtype;

    typedef typename SolverType::prec_raw_type    prec_raw_type;
    typedef typename SolverType::prec_type        prec_type;


    //! Constructor
    /*!
      \param dataType
      \param velocity FE space
      \param pressure FE space
      \param communicator
    */
    Oseen( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& uFESpace,
           FESpace<Mesh, EpetraMap>& pFESpace,
           boost::shared_ptr<Epetra_Comm>& comm,
           const int                 lagrangeMultiplier = 0);

    Oseen( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& uFESpace,
           FESpace<Mesh, EpetraMap>& pFESpace,
           boost::shared_ptr<Epetra_Comm>& comm,
           const EpetraMap           monolithicMap,
           const UInt                offset=0);

    /*!
      \param dataType
      \param lagrangeMultipliers (lagrange multipliers for the flux problem with rufaec flag)
      \param velocity FE space
      \param pressure FE space
      \param communicator
    */
    Oseen( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& uFESpace,
           FESpace<Mesh, EpetraMap>& pFESpace,
           std::vector<int> const&   lagrangeMultipliers,
           boost::shared_ptr<Epetra_Comm>& comm );


    //! virtual destructor

    virtual ~Oseen();

    //! Update convective term, bc treatment and solve the linearized ns system

    virtual void iterate( bchandler_raw_type& bch );

    virtual void setUp        ( const GetPot& dataFile );

    virtual void buildSystem();

    virtual void updateRHS(vector_type const& rhs)
    {
        M_rhsNoBC = rhs;
        M_rhsNoBC.GlobalAssemble();
    }

    virtual void updateSystem(const double       alpha,
                              const vector_type& betaVec,
                              const vector_type& sourceVec
                              );
    virtual void updateSystem(const double       alpha,
                              const vector_type& betaVec,
                              const vector_type& sourceVec,
                              matrix_ptrtype matrix
                              );

    void updateStab( matrix_type& matrFull );

//     void updateLinearSystem(double       alpha,
//                             vector_type& betaVec,
//                             vector_type& rhs
//                             );

    void initialize( const Function&, const Function&  );
    void initialize( const vector_type& u0, const vector_type& p0);

    void initialize( const vector_type& velAndPressure);

    //! returns the local solution vector

    const vector_ptrtype& solution() const {return M_sol;}

    //! returns the local residual vector

    const vector_type& residual() const {return M_residual;}

    //! reduce the local solution solution in global vectors

    void reduceSolution( Vector& u,
                         Vector& p );

    void reduceResidual( Vector& res);


    FESpace<Mesh, EpetraMap>& velFESpace()   {return M_uFESpace;}
    FESpace<Mesh, EpetraMap>& pressFESpace() {return M_pFESpace;}

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

    //! Returns the  Post Processing  stuff
    PostProc<Mesh>& post_proc();

    // Set up of post processing structures
    void post_proc_set_area();

    void post_proc_set_normal();

    void post_proc_set_phi();

    // compute area on a boundary face with given flag
    Real area(const EntityFlag& flag);

    // compute flux on a boundary face with given flag and a given solution
    Real flux(const EntityFlag& flag, const vector_type& sol);

    // compute average pressure on a boundary face with given flag and a given solution
    Real pressure(const EntityFlag& flag, const vector_type& sol);

    //! Get the lagrange multiplier related to a flux imposed on a given part of the boundary.
    /*!
     * @param Flag flag of the boundary face associated with the flux and the Lagrange multiplier we want.
     * @param BC BChandler containing the boundary conditions of the problem.
     * @return Lagrange multiplier
     */
    Real LagrangeMultiplier( const EntityFlag& Flag, bchandler_raw_type& BC );

    // compute flux on a boundary face with given flag
    Real flux(const EntityFlag& flag);

    // compute average pressure on a boundary face with given flag
    Real pressure(const EntityFlag& flag);

    //! Get the lagrange multiplier related to a flux imposed on a given part of the boundary.
    /*!
     * @param Flag flag of the boundary face associated with the flux and the Lagrange multiplier we want.
     * @param BC BChandler containing the boundary conditions of the problem.
     * @param Solution vector containing the solution of the problem (and also the Lagrange multipliers at the end).
     * @return Lagrange multiplier
     */
    Real LagrangeMultiplier( const EntityFlag& Flag, bchandler_raw_type& BC, const vector_type& Solution );

    // return the density and the viscosity of the fluid
    Real density()   const { return M_data.density(); }
    Real viscosity() const { return M_data.viscosity(); }

    //! Postprocessing
    void postProcess(bool _writeMesh = false);

    void resetPrec(bool reset = true) { if (reset) M_linearSolver.precReset(); }

    // as for now resetting stabilization matrix at the same time as the preconditioner
    void resetStab() { M_resetStab = true; }

    EpetraMap   const& getMap()       const { return M_localMap; }

    boost::shared_ptr<Epetra_Comm> const& comm()      const {return M_Displayer.comm();}

    Displayer   const& getDisplayer() const { return M_Displayer; }

    void recomputeMatrix(bool const recomp){M_recomputeMatrix = recomp;}

    matrix_type& matrNoBC()
        {
            return *M_matrNoBC;
        }
    matrix_type& matrMass()
        {
            return *M_matrMass;
        }

    bool          getIsDiagonalBlockPrec() {return M_isDiagonalBlockPrec;}
    void          setBlockPreconditioner(matrix_ptrtype blockPrec);
    void          getFluidMatrix( matrix_type& matrFull );
    void          updateUn( )                              {*M_un = *M_sol;}
    void          updateUn(const vector_type& sol )        {*M_un = sol;}// for the monolithic
    void          setupPostProc(const EntityFlag& flag, const Mesh meshPart);
protected:

    UInt dim_u() const           { return M_uFESpace.dim(); }
    UInt dim_p() const           { return M_pFESpace.dim(); }

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
    FESpace<Mesh, EpetraMap>&      M_uFESpace;
    FESpace<Mesh, EpetraMap>&      M_pFESpace;

    //! MPI communicator
    Displayer                      M_Displayer;

    EpetraMap                      M_localMap;

    //! mass matrix

    matrix_ptrtype                   M_matrMass;

    //! mass matrix

    matrix_ptrtype                   M_matrMassPr;

    //! Stokes matrix: mu*stiff
    matrix_ptrtype                   M_matrStokes;

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
    vector_ptrtype                    M_sol;

    //! residual

    vector_type                    M_residual;

    SolverType                     M_linearSolver;

    bool                           M_steady;

    //! Postprocessing class
    boost::shared_ptr<PostProc<Mesh> >                 M_post_proc;

    //! Stabilization
    bool                           M_stab;
    bool                           M_reuseStab;
    bool                           M_resetStab;
    int                            M_iterReuseStab;

    details::
    IPStabilization<Mesh, Dof>     M_ipStab;
    Real                           M_gammaBeta;
    Real                           M_gammaDiv;
    Real                           M_gammaPress;

    const Function*                M_betaFct;

    bool                           M_divBetaUv;

    bool                           M_stiffStrain;

    //
    double                         M_diagonalize;

    UInt                           M_count;

    bool                           M_recomputeMatrix;

    bool                           M_isDiagonalBlockPrec;

  //private:

    //! Elementary matrices and vectors
    ElemMat                        M_elmatStiff;      // velocity Stokes
    ElemMat                        M_elmatMass;       // velocity mass
    ElemMat                        M_elmatP;          // (p,q) bloc for preconditioners
    ElemMat                        M_elmatDiv;
    ElemMat                        M_elmatGrad;
    ElemVec                        M_elvec;           // Elementary right hand side
    matrix_ptrtype                 M_blockPrec;
    ElemVec                        M_wLoc;
    ElemVec                        M_uLoc;
    boost::shared_ptr<vector_type> M_un;
}; // class Oseen



//
// IMPLEMENTATION
//


template<typename Mesh, typename SolverType>
Oseen<Mesh, SolverType>::
Oseen( const data_type&          dataType,
       FESpace<Mesh, EpetraMap>& uFESpace,
       FESpace<Mesh, EpetraMap>& pFESpace,
       boost::shared_ptr<Epetra_Comm>&              comm,
       const int                 lagrangeMultiplier):
    M_data                   ( dataType ),
    M_uFESpace               ( uFESpace ),
    M_pFESpace               ( pFESpace ),
    M_Displayer              ( comm ),
    M_localMap               ( M_uFESpace.map() + M_pFESpace.map() + lagrangeMultiplier),
    M_matrMass               ( ),
    M_matrMassPr             ( ),
    M_matrStokes             ( ),
//    M_matrFull               ( ),
    M_matrNoBC               ( ),
    M_matrStab               ( ),
    M_rhsNoBC                ( M_localMap ),
    M_rhsFull                ( M_localMap ),
    M_sol                    ( new vector_type(M_localMap) ),
    M_residual               ( M_localMap ),
    M_linearSolver           ( comm ),
    M_post_proc              ( new PostProc<Mesh>(M_uFESpace.mesh(),
                                                  &M_uFESpace.feBd(), &M_uFESpace.dof(),
                                                  &M_pFESpace.feBd(), &M_pFESpace.dof(), M_localMap )),
    M_stab                   ( false ),
    M_reuseStab              ( false ),
    M_resetStab              ( false ),
    M_iterReuseStab          ( -1 ),
    M_ipStab                 ( M_uFESpace.mesh(),
                               M_uFESpace.dof(), M_uFESpace.refFE(),
                               M_uFESpace.feBd(), M_uFESpace.qr(),
                               0., 0., 0.,
                               M_data.viscosity() ),
    M_betaFct                ( 0 ),
    M_divBetaUv              ( false ),
    M_stiffStrain            ( false ),
    M_diagonalize            ( false ),
    M_count                  ( 0 ),
    //    M_updated                ( false ),
    M_recomputeMatrix        ( false ),
    M_elmatStiff             ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatMass              ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatP                 ( M_pFESpace.fe().nbNode, 1, 1 ),
    M_elmatDiv               ( M_pFESpace.fe().nbNode, 1, 0, M_uFESpace.fe().nbNode, 0, nDimensions ),
    M_elmatGrad              ( M_uFESpace.fe().nbNode, nDimensions, 0, M_pFESpace.fe().nbNode, 0, 1 ),
    M_elvec                  ( M_uFESpace.fe().nbNode, nDimensions ),
    M_blockPrec              (),
    M_wLoc                   ( M_uFESpace.fe().nbNode, nDimensions ),
    M_uLoc                   ( M_uFESpace.fe().nbNode, nDimensions ),
    M_un                     (new vector_type(M_localMap))
{
    M_stab = (&M_uFESpace.refFE() == &M_pFESpace.refFE());
}

template<typename Mesh, typename SolverType>
Oseen<Mesh, SolverType>::
Oseen( const data_type&          dataType,
       FESpace<Mesh, EpetraMap>& uFESpace,
       FESpace<Mesh, EpetraMap>& pFESpace,
       boost::shared_ptr<Epetra_Comm>&              comm ,
       EpetraMap                 monolithicMap,
       UInt                      /*offset*/):
    M_data                   ( dataType ),
    M_uFESpace               ( uFESpace ),
    M_pFESpace               ( pFESpace ),
    M_Displayer              ( comm ),
    M_localMap               ( monolithicMap ),
    M_matrMass               ( ),
    M_matrStokes             ( ),
    //    M_matrFull               ( ),
    M_matrNoBC               ( ),
    M_matrStab               ( ),
    M_rhsNoBC                ( M_localMap ),
    M_rhsFull                ( M_localMap ),
    M_sol                    ( ),
    M_residual               ( M_localMap ),
    M_linearSolver           ( comm ),
    M_post_proc              ( new PostProc<Mesh>(M_uFESpace.mesh(),
                                                  &M_uFESpace.feBd(), &M_uFESpace.dof(),
                                                  &M_pFESpace.feBd(), &M_pFESpace.dof(), M_localMap )),
    M_stab                   ( false ),
    M_reuseStab              ( false ),
    M_resetStab              ( false ),
    M_iterReuseStab          ( -1 ),
    M_ipStab                 ( M_uFESpace.mesh(),
                               M_uFESpace.dof(), M_uFESpace.refFE(),
                               M_uFESpace.feBd(), M_uFESpace.qr(),
                               0., 0., 0.,
                               M_data.viscosity() ),
    M_betaFct                ( 0 ),
    M_divBetaUv              ( false ),
    M_stiffStrain            ( false ),
    M_diagonalize            ( false ),
    M_count                  ( 0 ),
    //    M_updated                ( false ),
    M_recomputeMatrix        ( false ),
    M_elmatStiff             ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatMass              ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatP                 ( M_pFESpace.fe().nbNode, 1, 1 ),
    M_elmatDiv               ( M_pFESpace.fe().nbNode, 1, 0, M_uFESpace.fe().nbNode, 0, nDimensions ),
    M_elmatGrad              ( M_uFESpace.fe().nbNode, nDimensions, 0, M_pFESpace.fe().nbNode, 0, 1 ),
    M_elvec                  ( M_uFESpace.fe().nbNode, nDimensions ),
    M_blockPrec              (),
    M_wLoc                   ( M_uFESpace.fe().nbNode, nDimensions ),
    M_uLoc                   ( M_uFESpace.fe().nbNode, nDimensions ),
    M_un                     (new vector_type(M_localMap))
{
    M_stab = (&M_uFESpace.refFE() == &M_pFESpace.refFE());
}

template<typename Mesh, typename SolverType>
Oseen<Mesh, SolverType>::
Oseen( const data_type&          dataType,
       FESpace<Mesh, EpetraMap>& uFESpace,
       FESpace<Mesh, EpetraMap>& pFESpace,
       std::vector<int> const&   lagrangeMultipliers,
       boost::shared_ptr<Epetra_Comm>& comm ):
    M_data                   ( dataType ),
    M_uFESpace               ( uFESpace ),
    M_pFESpace               ( pFESpace ),
    M_Displayer              ( comm ),
    M_localMap               ( M_uFESpace.map() + M_pFESpace.map() + lagrangeMultipliers ),
    M_matrMass               ( ),
    M_matrStokes             ( ),
    //    M_matrFull               ( ),
    M_matrNoBC               ( ),
    M_matrStab               ( ),
    M_rhsNoBC                ( M_localMap ),
    M_rhsFull                ( M_localMap ),
    M_sol                    ( new vector_type(M_localMap) ),
    M_residual               ( M_localMap ),
    M_linearSolver           ( ),
    M_post_proc              ( new PostProc<Mesh>(M_uFESpace.mesh(),
                                                  &M_uFESpace.feBd(), &M_uFESpace.dof(),
                                                  &M_pFESpace.feBd(), &M_pFESpace.dof(), M_localMap )),
    M_stab                   ( false ),
    M_reuseStab              ( false ),
    M_resetStab              ( false ),
    M_iterReuseStab          ( -1 ),
    M_ipStab                 ( M_uFESpace.mesh(),
                               M_uFESpace.dof(), M_uFESpace.refFE(),
                               M_uFESpace.feBd(), M_uFESpace.qr(),
                               0., 0., 0.,
                               M_data.viscosity() ),
    M_betaFct                ( 0 ),
    M_divBetaUv              ( false ),
    M_stiffStrain            ( false ),
    M_diagonalize            ( false ),
    M_count                  ( 0 ),
    //    M_updated                ( false ),
    M_recomputeMatrix        ( false ),
    M_elmatStiff             ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatMass              ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatP                 ( M_pFESpace.fe().nbNode, 1, 1 ),
    M_elmatDiv               ( M_pFESpace.fe().nbNode, 1, 0, M_uFESpace.fe().nbNode, 0, nDimensions ),
    M_elmatGrad              ( M_uFESpace.fe().nbNode, nDimensions, 0, M_pFESpace.fe().nbNode, 0, 1 ),
    M_elvec                  ( M_uFESpace.fe().nbNode, nDimensions ),
    M_blockPrec              (),
    M_wLoc                   ( M_uFESpace.fe().nbNode, nDimensions ),
    M_uLoc                   ( M_uFESpace.fe().nbNode, nDimensions ),
    M_un                     (new vector_type(M_localMap))
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

    M_linearSolver.setUpPrec(dataFile, "fluid/prec");
    M_linearSolver.setDataFromGetPot( dataFile, "fluid/solver" );

    M_steady        = dataFile( "fluid/miscellaneous/steady",        0  );

    M_gammaBeta     = dataFile( "fluid/ipstab/gammaBeta",            0. );
    M_gammaDiv      = dataFile( "fluid/ipstab/gammaDiv",             0. );
    M_gammaPress    = dataFile( "fluid/ipstab/gammaPress",           0. );
    M_reuseStab     = dataFile( "fluid/ipstab/reuse",               false);
    M_iterReuseStab = dataFile( "fluid/ipstab/max_iter_reuse", static_cast<int> ( M_linearSolver.MaxIter() * 8./10. ) );

    M_divBetaUv   = dataFile( "fluid/space_discretization/div_beta_u_v",false); // Energetic stabilization term
    M_stiffStrain = dataFile( "fluid/space_discretization/stiff_strain",false); // Enable grad( u )^T in stress tensor
    M_diagonalize = dataFile( "fluid/space_discretization/diagonalize",  1. );
    M_isDiagonalBlockPrec = dataFile( "fluid/diagonalBlockPrec",  false );


    //    M_linearSolver.setAztecooPreconditioner( dataFile, "fluid/solver" );

    M_ipStab.setGammaBeta (M_gammaBeta);
    M_ipStab.setGammaDiv  (M_gammaDiv);
    M_ipStab.setGammaPress(M_gammaPress);

}




template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::buildSystem()
{


    M_matrMass.reset  ( new matrix_type(M_localMap) );
    M_matrStokes.reset( new matrix_type(M_localMap) );

    M_Displayer.leaderPrint("  F-  Computing constant matrices ...          ");

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

    if(M_isDiagonalBlockPrec == true)
        {
            M_blockPrec.reset(new matrix_type(M_localMap));
        }
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

        // stiffness
        chronoStiff.start();
        if ( M_stiffStrain )
            stiff_strain( 2.0*M_data.viscosity(), M_elmatStiff, M_uFESpace.fe() );
        else
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
            // stiffness
            chronoStiffAssemble.start();
            if ( M_isDiagonalBlockPrec == true )
            {
                assembleMatrix( *M_blockPrec,
                                M_elmatStiff,
                                M_uFESpace.fe(),
                                M_uFESpace.fe(),
                                M_uFESpace.dof(),
                                M_uFESpace.dof(),
                                iComp, iComp,
                                iComp*velTotalDof, iComp*velTotalDof);
            }
            else
            {
                if ( M_stiffStrain ) // sigma = 0.5 * mu (grad( u ) + grad ( u )^T)
                {
                    for ( UInt jComp = 0; jComp < nbCompU; jComp++ )
                    {
                        assembleMatrix( *M_matrStokes,
                                        M_elmatStiff,
                                        M_uFESpace.fe(),
                                        M_uFESpace.fe(),
                                        M_uFESpace.dof(),
                                        M_uFESpace.dof(),
                                        iComp, jComp,
                                        iComp*velTotalDof, jComp*velTotalDof);

                    }
                }
                else // sigma = mu grad( u )
                {
                    assembleMatrix( *M_matrStokes,
                                    M_elmatStiff,
                                    M_uFESpace.fe(),
                                    M_uFESpace.fe(),
                                    M_uFESpace.dof(),
                                    M_uFESpace.dof(),
                                    iComp, iComp,
                                    iComp*velTotalDof, iComp*velTotalDof);
                }
            }
            chronoStiffAssemble.stop();

            // mass
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





    //    for (UInt ii = nDimensions*dim_u(); ii < nDimensions*dim_u() + dim_p(); ++ii)
    //  M_matrStokes->set_mat_inc( ii ,ii, 0. ); not scalable!!!

    if(M_isDiagonalBlockPrec == true)
        {
            M_blockPrec->GlobalAssemble();
            *M_matrStokes += *M_blockPrec;
        }
    comm()->Barrier();

    chrono.stop();
    M_Displayer.leaderPrintMax( "done in " , chrono.diff());


    M_Displayer.leaderPrint( "  F-  Finalizing the matrices ...              ");

    chrono.start();

    M_matrStokes->GlobalAssemble();
    M_matrMass->GlobalAssemble();

    chrono.stop();
    M_Displayer.leaderPrintMax("done in " , chrono.diff() );

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
void Oseen<Mesh, SolverType>::
initialize( const Function& u0, const Function& p0 )
{
     vector_type u(M_uFESpace.map());
     M_uFESpace.interpolate(u0, u, M_data.dataTime()->getTime());

     vector_type p(M_pFESpace.map());
     M_pFESpace.interpolate(p0, p, M_data.dataTime()->getTime());

     initialize(u, p);
}


template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::
initialize( const vector_type& u0, const vector_type& p0 )
{

    *M_sol = u0;
    *M_un = u0;
    M_sol->add(p0, nDimensions*M_uFESpace.dof().numTotalDof());
    M_un->add(p0, nDimensions*M_uFESpace.dof().numTotalDof());

}

template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::
initialize( const vector_type& velAndPressure)
{

    *M_sol = velAndPressure;
    *M_un = velAndPressure;

}



template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::
updateSystem(const double       alpha,
             const vector_type& betaVec,
             const vector_type& sourceVec
             )
{
    if (M_matrNoBC.get())
        M_matrNoBC.reset(new matrix_type(M_localMap, M_matrNoBC->getMeanNumEntries() ));
    else
        M_matrNoBC.reset(new matrix_type(M_localMap));
    updateSystem( alpha, betaVec, sourceVec, M_matrNoBC);
    if(alpha != 0.)
        M_matrNoBC->GlobalAssemble();

}

template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::
updateSystem(const double       alpha,
             const vector_type& betaVec,
             const vector_type& sourceVec,
             matrix_ptrtype matrNoBC)
{
    Chrono chrono;

    // clearing pressure mass matrix in case we need it in removeMean;
    M_matrMassPr.reset( );


    M_Displayer.leaderPrint("  F-  Updating mass term on right hand side... ");

    chrono.start();

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();
//    UInt pressTotalDof = M_pFESpace.dof().numTotalDof();

    // Right hand side for the velocity at time

    updateRHS(sourceVec);

    chrono.stop();

    M_Displayer.leaderPrintMax("done in ", chrono.diff());


    //    M_updated = false;

//

    if (M_recomputeMatrix)
        buildSystem();

    M_Displayer.leaderPrint( "  F-  Copying the matrices ...                 ");

    chrono.start();

    if(M_isDiagonalBlockPrec == true)
        {
            matrix_ptrtype tmp(M_blockPrec);
            M_blockPrec.reset(new matrix_type(M_localMap, M_blockPrec->getMeanNumEntries() ));
            *M_blockPrec += *tmp;
        }


    chrono.stop();
    M_Displayer.leaderPrintMax( "done in " , chrono.diff() );


    UInt nbCompU       = nDimensions;

    //! managing the convective term

    double normInf;
    betaVec.NormInf(&normInf);

    if (normInf != 0.)
    {
        M_Displayer.leaderPrint("  F-  Sharing convective term ...              ");
        chrono.start();

        // vector with repeated nodes over the processors

        vector_type betaVecRep(betaVec, Repeated );
        vector_type unRep(*M_un, Repeated);

        chrono.stop();

        M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
        M_Displayer.leaderPrint("  F-  Updating the convective terms ...        ");
        chrono.start();

        for ( UInt iVol = 1; iVol<= M_uFESpace.mesh()->numVolumes(); ++iVol )
        {

            M_pFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) ); // just to provide the id number in the assem_mat_mixed
            M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) ); //as updateFirstDer

            M_elmatStiff.zero();

            UInt eleID = M_uFESpace.fe().currentLocalId();
            // Non linear term, Semi-implicit approach
            // M_elvec contains the velocity values in the nodes
            for ( UInt iNode = 0 ; iNode < M_uFESpace.fe().nbNode ; iNode++ )
            {
                UInt  iloc = M_uFESpace.fe().patternFirst( iNode );
                for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
                {
                    UInt ig = M_uFESpace.dof().localToGlobal( eleID, iloc + 1 ) + iComp*dim_u();
                    M_elvec.vec()[ iloc + iComp*M_uFESpace.fe().nbNode ] = betaVecRep[ig]; // BASEINDEX + 1

                    M_uLoc.vec() [ iloc + iComp * M_uFESpace.fe().nbNode ] = (unRep)(ig);
                    M_wLoc.vec() [ iloc + iComp * M_uFESpace.fe().nbNode ] = (unRep)(ig)-betaVecRep(ig);
                }
            }


            // ALE term: - rho div w u v
            mass_divw( -M_data.density(), M_wLoc,  M_elmatStiff , M_uFESpace.fe(), 0, 0, nbCompU );//ADDED

            // ALE stab implicit: 0.5 rho div u w v
            mass_divw( 0.5*M_data.density(), M_uLoc,  M_elmatStiff , M_uFESpace.fe(), 0, 0, nbCompU );//ADDED


            // Stabilising term: div u^n u v
            if ( M_divBetaUv )
                mass_divw( 0.5*M_data.density(), M_elvec, M_elmatStiff, M_uFESpace.fe(), 0, 0, nbCompU );

            // compute local convective terms
            advection( M_data.density(), M_elvec, M_elmatStiff, M_uFESpace.fe(), 0, 0, nbCompU );

            // loop on components
            for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
            {
                // compute local convective term and assembling
                // grad( 0, M_elvec, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );
                // grad( 1, M_elvec, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );
                // grad( 2, M_elvec, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );

                assembleMatrix( *matrNoBC,
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
        M_Displayer.leaderPrintMax( "done in " , chrono.diff() );

        if ( M_stab && (M_resetStab || !M_reuseStab || (M_matrStab.get() == 0) ) )
        {
            M_Displayer.leaderPrint("  F-  Updating the stabilization terms ...     ");
            chrono.start();
            M_matrStab.reset  ( new matrix_type(M_localMap) );
            M_ipStab.apply( *M_matrStab, betaVecRep, false );
            M_matrStab->GlobalAssemble();
            M_resetStab = false;
            chrono.stop();
            M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
        }

    }
    else
        {
            if (M_stab)
                {
                    M_Displayer.leaderPrint("  F-  Updating the stabilization terms ...     ");
                    chrono.start();

                    if ( M_resetStab || !M_reuseStab || (M_matrStab.get() == 0) )
                        {
                            M_matrStab.reset  ( new matrix_type(M_localMap) );
                            M_ipStab.apply( *M_matrStab, betaVec, false );
                            M_matrStab->GlobalAssemble();
                            M_resetStab = false;

                        }
                    else
                        {
                            M_Displayer.leaderPrint("reusing stab. ");
                        }
                    chrono.stop();
                    M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
                }
        }

    if (alpha != 0. )
        {
            *matrNoBC += *M_matrMass*alpha;
            if(M_isDiagonalBlockPrec == true)
                {
                    matrNoBC->GlobalAssemble();
                    *M_blockPrec += *matrNoBC;
                    matrix_type tmp(*matrNoBC);
                    matrNoBC.reset(new matrix_type(M_localMap, tmp.getMeanNumEntries()));
                    *matrNoBC += tmp;
                    M_blockPrec->GlobalAssemble();
                }
        }
    *matrNoBC += *M_matrStokes;


}


template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::setBlockPreconditioner( matrix_ptrtype blockPrec )
{
    //    blockPrec.reset(new matrix_type(M_monolithicMap, M_solid->getMatrixPtr()->getMeanNumEntries()));
    *blockPrec += *M_blockPrec;
}

template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::updateStab( matrix_type& matrFull )
{

    if (M_stab)
    {
        matrFull += *M_matrStab;
    }

}

template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::getFluidMatrix( matrix_type& matrFull )
{
    M_matrNoBC->GlobalAssemble();
    matrFull += *M_matrNoBC;
}

template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::iterate( bchandler_raw_type& bch )
{

    Chrono chrono;

    // matrix and vector assembling communication
    M_Displayer.leaderPrint("  F-  Updating the boundary conditions ...     ");

    chrono.start();


    M_matrNoBC->GlobalAssemble();

    matrix_ptrtype matrFull( new matrix_type( M_localMap, M_matrNoBC->getMeanNumEntries()));

    updateStab(*matrFull);
    getFluidMatrix(*matrFull);

    vector_type rhsFull (M_rhsNoBC);

//     matrFull.reset(new matrix_type(*M_matrNoBC));
//     M_rhsFull = M_rhsNoBC;

    chrono.stop();

    M_Displayer.leaderPrintMax("done in ", chrono.diff() );

    // boundary conditions update
    M_Displayer.leaderPrint("  F-  Applying boundary conditions ...         ");

    chrono.start();
    applyBoundaryConditions( *matrFull, rhsFull, bch);

    matrFull->GlobalAssemble();
    chrono.stop();

    M_Displayer.leaderPrintMax("done in " , chrono.diff());

    //matrFull->spy("matrFull");
    // solving the system
    M_linearSolver.setMatrix(*matrFull);

    int numIter = M_linearSolver.solveSystem( rhsFull, *M_sol, matrFull );

    // if the preconditioner has been rese the stab terms are to be updated
    if( numIter < 0 || numIter > M_iterReuseStab )
    {
        resetStab();
    }

    M_residual  = M_rhsNoBC;
    M_residual -= *M_matrNoBC*(*M_sol);

    //M_residual.spy("residual");
} // iterate()



template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::reduceSolution( Vector& u,
                                              Vector& p )
{
    vector_type vel(*M_sol, 0);

    if (false/*S_verbose*/)
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

    if (false/*S_verbose*/)
    {
        for ( UInt iDof = 0; iDof < nDimensions*dim_u(); ++iDof )
        {
            res[ iDof ] = vel[ iDof + 1 ]; // BASEINDEX + 1
        }

    }
}


template<typename Mesh, typename SolverType>
Real Oseen<Mesh, SolverType>::removeMean( vector_type& x )
{

    Chrono chrono;
    chrono.start();

    UInt const nbCompU (nDimensions);
    UInt const velTotalDof (M_uFESpace.dof().numTotalDof());


    if ( M_matrMassPr.get() == 0)
        M_matrMassPr.reset  ( new matrix_type(M_localMap) );

    for ( UInt iVol = 1; iVol <= M_uFESpace.mesh()->numVolumes(); iVol++ )
    {
        chrono.start();
        M_pFESpace.fe().update( M_pFESpace.mesh()->volumeList( iVol ) ); // just to provide the id number in the assem_mat_mixed



        M_elmatP.zero();
        // mass
        chrono.start();
        mass( 1, M_elmatP, M_pFESpace.fe(), 0, 0, nDimensions );
        chrono.stop();

        chrono.start();
        assembleMatrix( *M_matrMassPr,
                        M_elmatP,
                        M_pFESpace.fe(),
                        M_pFESpace.fe(),
                        M_pFESpace.dof(),
                        M_pFESpace.dof(),
                        nbCompU, nbCompU,
                        nbCompU*velTotalDof, nbCompU*velTotalDof);
        chrono.stop();
    }

    M_matrMassPr->GlobalAssemble();

    vector_type ones(*M_sol);
    ones = 1.0;

    Real mean;
    mean = ones* ( M_matrMassPr * x);
    x += (-mean);

    ASSERT( abs(ones* ( M_matrMassPr * x)) < 1e-9 , "after removeMean the mean pressure should be zero!");

    return mean;


} // removeMean()

template<typename Mesh, typename SolverType>
void Oseen<Mesh, SolverType>::applyBoundaryConditions( matrix_type&        matrix,
                                                       vector_type&        rhs,
                                                       bchandler_raw_type& BCh )
{

    // M_rhsFull = M_rhsNoBC;

    // BC manage for the velocity
    if ( !BCh.bdUpdateDone() ||  M_recomputeMatrix )
    {
        BCh.bdUpdate( *M_uFESpace.mesh(), M_uFESpace.feBd(), M_uFESpace.dof() );
    }

    //vector_type rhsFull(rhs, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately
    vector_type rhsFull(rhs, Unique); // ignoring non-local entries, Otherwise they are summed up lately

    bcManage( matrix, rhsFull, *M_uFESpace.mesh(), M_uFESpace.dof(), BCh, M_uFESpace.feBd(), 1.,
              M_data.dataTime()->getTime() );

    rhs = rhsFull;

    if ( BCh.hasOnlyEssential() && M_diagonalize )
    {
        matrix.diagonalize( nDimensions*dim_u(),
                            M_diagonalize,
                            rhs,
                            0.);
    }

} // applyBoundaryCondition



template<typename Mesh, typename SolverType> void
Oseen<Mesh, SolverType>::setupPostProc(const EntityFlag& flag, const Mesh meshPart)
{
    M_post_proc.reset(new PostProc<Mesh>(M_uFESpace.mesh(), &M_uFESpace.feBd(), &M_uFESpace.dof(), &M_pFESpace.feBd(), &M_pFESpace.dof(), M_localMap ));
}


// Returns the Post Processing structure
template<typename Mesh, typename SolverType>
PostProc<Mesh>&
Oseen<Mesh, SolverType>::post_proc()
{
    return *M_post_proc;
}

// Set up of post processing structures
template<typename Mesh, typename SolverType>
void
Oseen<Mesh, SolverType>::post_proc_set_area()
{
  M_post_proc->set_area();
}

template<typename Mesh, typename SolverType>
void
Oseen<Mesh, SolverType>::post_proc_set_normal()
{
  M_post_proc->set_normal();
}

template<typename Mesh, typename SolverType>
void
Oseen<Mesh, SolverType>::post_proc_set_phi()
{
    M_post_proc->set_phi();
}



//! Computes the flux on a given part of the boundary
template<typename Mesh, typename SolverType> Real
Oseen<Mesh, SolverType>::flux(const EntityFlag& flag)
{
    return flux(flag, *M_sol);
}

//! Computes the pressure on a given part of the boundary
template<typename Mesh, typename SolverType> Real
Oseen<Mesh, SolverType>::pressure(const EntityFlag& flag)
{
  return pressure(flag, *M_sol);
}

template<typename Mesh, typename SolverType> Real
Oseen<Mesh, SolverType>::LagrangeMultiplier( const EntityFlag& Flag, bchandler_raw_type& BC )
{
  return LagrangeMultiplier( Flag, BC, *M_sol );
}

//! Computes the flux on a given part of the boundary
template<typename Mesh, typename SolverType> Real
Oseen<Mesh, SolverType>::flux(const EntityFlag& flag, const vector_type& sol)
{
  vector_type velAndPressure(sol, Repeated);
  vector_type vel(this->M_uFESpace.map(), Repeated);
  vel.subset(velAndPressure);

  return M_post_proc->flux(vel, flag);
}

//! Computes the pressure on a given part of the boundary
template<typename Mesh, typename SolverType> Real
Oseen<Mesh, SolverType>::pressure(const EntityFlag& flag, const vector_type& sol)
{
  vector_type velAndPressure(sol, Repeated);
  vector_type press(this->M_pFESpace.map(), Repeated);
  press.subset(velAndPressure, this->M_uFESpace.dim()*this->M_uFESpace.fieldDim());

  // third argument is 1, to use the pressure finite element space (see PostProc docs)
  return M_post_proc->average(press, flag, 1)[0];
}

template<typename Mesh, typename SolverType> Real
Oseen<Mesh, SolverType>::LagrangeMultiplier( const EntityFlag&         Flag,
                                                   bchandler_raw_type& BC,
                                             const vector_type&        Solution )
{
    // Create a list of Flux BCName
    std::vector< BCName > FluxBCVector = BC.getBCWithType( Flux );
    BCName FluxBCName = BC.GetBCWithFlag( Flag ).name();

    // Create a Repeated vector for looking to the lambda
    vector_type VelocityPressureLambda(Solution, Repeated);

    // Find the index associated to the correct Lagrange multiplier
    for ( UInt lmIndex(0); lmIndex < static_cast <UInt> ( FluxBCVector.size() ); ++lmIndex )
        if ( FluxBCName.compare( FluxBCVector[lmIndex] ) == 0 )
            return VelocityPressureLambda[3 * M_uFESpace.dof().numTotalDof() + M_pFESpace.dof().numTotalDof() + 1 + lmIndex];

    // If lmIndex has not been found a warning message is printed
    std::cout << "!!! Warning - Lagrange multiplier for Flux BC not found!" << std::endl;
    return 0;
}

//! Computes the area on a given part of the boundary
template<typename Mesh, typename SolverType> Real
Oseen<Mesh, SolverType>::area(const EntityFlag& flag) {

  return M_post_proc->area(flag);
}

// Postprocessing
template <typename Mesh, typename SolverType>
void
Oseen<Mesh, SolverType>::postProcess(bool /*_writeMesh*/)
{
    std::ostringstream index;
    std::ostringstream indexMe;
    std::string name;
    std::string me;


    indexMe << M_Displayer.comm()->MyPID();

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

    vector_type velAndPressure(*M_sol, Repeated);
    vector_type res(M_residual, Repeated);


//         if ( fmod( float( M_count ), float( M_data.verbose() ) ) == 0.0 )
//         {
    M_Displayer.leaderPrint( "  F-  Post-processing " );

    index << std::setfill('0') << std::setw(3);
    index << ( M_count / M_data.verbose() );
    name = index.str();

    PhysVectUnknown<Vector> u(nDimensions*dim_u());
    ScalUnknown<Vector>     p(dim_p());

    reduceSolution(u, p);

//     if (M_me == 0)
//     {
        // postprocess data file for medit
//         wr_medit_ascii_scalar( "vel_x." + name + ".bb", u.giveVec(),
//                                M_data.mesh()->numGlobalVertices() );
//         wr_medit_ascii_scalar( "vel_y." + name + ".bb", u.giveVec() + this->dim_u(),
//                                M_data.mesh()->numGlobalVertices() );
//         wr_medit_ascii_scalar( "vel_z." + name + ".bb", u.giveVec()+2*this->dim_u(),
//                                M_data.mesh()->numGlobalVertices() );
//         wr_medit_ascii_scalar( "press." + name + ".bb", p.giveVec(),
//                                p.size() );

    	// double dt = M_data.getTimeStep();


       writeMesh("vel_x." + me + "." + name + ".mesh", *M_uFESpace.mesh());
       writeMesh("resf_x." + me +  "." + name + ".mesh", *M_uFESpace.mesh());

       writeMesh("vel_y." + me + "." + name + ".mesh", *M_uFESpace.mesh());
       writeMesh("resf_y." + me +  "." + name + ".mesh", *M_uFESpace.mesh());

       writeMesh("vel_z." + me + "." + name + ".mesh", *M_uFESpace.mesh());
       writeMesh("resf_z." + me +  "." + name + ".mesh", *M_uFESpace.mesh());

       writeMesh("press." + me +  "." + name + ".mesh", *M_uFESpace.mesh());


        meditSolutionWriter( "vel_x." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), velAndPressure, M_uFESpace.dof().numTotalDof()*0);
        meditSolutionWriter( "vel_y." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), velAndPressure, M_uFESpace.dof().numTotalDof()*1);
        meditSolutionWriter( "vel_z." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), velAndPressure, M_uFESpace.dof().numTotalDof()*2);

//         wr_medit_ascii2("vel_x." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );
//         wr_medit_ascii2("vel_y." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );
//         wr_medit_ascii2("vel_z." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );


//         wr_medit_ascii2("resf_x." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );
//         wr_medit_ascii2("resf_y." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );
//         wr_medit_ascii2("resf_z." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );

        meditSolutionWriter( "resf_x." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), res, M_uFESpace.dof().numTotalDof()*0);
        meditSolutionWriter( "resf_y." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), res, M_uFESpace.dof().numTotalDof()*1);
        meditSolutionWriter( "resf_z." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), res, M_uFESpace.dof().numTotalDof()*2);


        meditSolutionWriter( "press." + me + "." + name + ".bb",
                             *M_pFESpace.mesh(), velAndPressure, M_uFESpace.dof().numTotalDof()*3);


//    cout << ( "should do ln -s -f " + M_data.meshDir() + M_data.meshFile() +
//             " press." + name + "." + me + ".mesh" ).data()  ;

//         system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
//                   " press." + name + ".mesh" ).data() );
//         system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
//                   " vel_x." + name + ".mesh" ).data() );
//         system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
//                   " vel_y." + name + ".mesh" ).data() );
//         system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
//                   " vel_z." + name + ".mesh" ).data() );

//     }
//    }

        M_count++;

}



} // namespace LifeV


#endif //_OSEEN_H_
