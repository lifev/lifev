/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
  \file VenantKirchhofSolver.h
  \author M.A. Fernandez
  \date 6/2003
  \version 1.0

  \brief
  This file contains solvers for St. Venant-Kirchhof materials (linear for the moment)

*/
#ifndef _VENANTKIRCHHOFSOLVER_H_
#define _VENANTKIRCHHOFSOLVER_H_

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifealg/dataNewton.hpp>
#include <life/lifealg/newton.hpp>
//
#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>

//
#include "life/lifefem/FESpace.hpp"
#include <life/lifesolver/dataElasticStructure.hpp>



namespace LifeV
{


/*!
  \class VenantKirchhofSolver
  \brief
  This class solves the linear elastodynamics equations for a (only linear right now)
  St. Venant-Kirchoff material

*/
template <typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >

class VenantKirchhofSolver
{
public:

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> source_type;

    typedef BCHandler                             bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type> bchandler_type;


    typedef SolverType                    solver_type;

    typedef typename solver_type::matrix_type      matrix_type;
    typedef boost::shared_ptr<matrix_type>         matrix_ptrtype;
    typedef typename solver_type::vector_type      vector_type;

    typedef typename SolverType::prec_raw_type    prec_raw_type;
    typedef typename SolverType::prec_type        prec_type;

    typedef DataElasticStructure<Mesh>    data_type;



    //! Constructors
    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature rule
      \param BCh boundary conditions for the displacement
    */

    VenantKirchhofSolver( const data_type& data,
                          FESpace<Mesh, EpetraMap>&   FESpace,
                          BCHandler&       BCh,
                          Epetra_Comm&     comm,
                          UInt             offset=0);

    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature rule
    */

    VenantKirchhofSolver( const data_type& data,
                          FESpace<Mesh, EpetraMap>&   FESpace,
                          Epetra_Comm&     comm,
                          UInt       offset=0);

    VenantKirchhofSolver( const data_type& data,
                          FESpace<Mesh, EpetraMap>&   dFESpace,
                          Epetra_Comm&     comm,
                          EpetraMap&       monolithicMap,
                          UInt       offset=0,
                          boost::shared_ptr<FESpace<Mesh, EpetraMap> >   uFESpace=0);


    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
//    void updateSystem( source_type const& source );
    void updateSystem(vector_type& rhsFluidCoupling); // used for monolithic
    void updateSystem();

    void buildSystem(matrix_type & bigMatrixStokes); // used for monolithic
    void buildSystem();

    //! Solve the non-linear system
    void iterateMonolithic(vector_type& rhs, vector_type& sol, matrix_ptrtype prec); // used for monolithic
    void iterate( vector_type& sol );
    void iterate( bchandler_raw_type& bch );
    void iterateLin( bchandler_raw_type& bch );

    //! Output
    //void showMe( std::ostream& c = std::cout ) const;

    //! getters

    //! BCHandler getter and setter

//    LIFEV_DEPRECATED BCHandler const & BC_solid() const {return BCh_solid();}

    FESpace<Mesh, EpetraMap>& dFESpace(){return M_FESpace;}


    BCHandler const & BChandler() const {return M_BCh;}

    void setBC(BCHandler& BCd)   {M_BCh = &BCd;}

    //! residual getter
    vector_type& residual()             {return M_residual_d;}

    //! recur setter

    void setRecur(UInt recur) {_recur = recur;}

    void updateJacobian( vector_type& sol, int iter );

    //! solves the tangent problem for newton iterations
    void solveJac( vector_type&       step,
                   const vector_type& res,
                   double&            linear_rel_tol);
    //    void solveJac( const Vector& res, double& linear_rel_tol, Vector &step);
    //! solves the tangent problem with custom BC

    void solveJac( vector_type&       step,
                   const vector_type& res,
                   double&            linear_rel_tol,
                   bchandler_type&    BCd );

    void solveJacobian( const Real );
//    void solveJacobian( );
    void solveJacobian( const           Real,
                        bchandler_type& BCd);

//! evaluates residual for newton interations
    void evalResidual( vector_type &res, const vector_type& sol, int iter);

    void evalConstraintTensor();

    void setSourceTerm( source_type const& __s ) { M_source = __s; }

    source_type const& sourceTerm() const { return M_source; }

    vector_type& disp()        { return M_disp; }
    vector_type& vel()         { return M_vel; }
    vector_type& rhsWithoutBC() { return M_rhsNoBC; }

    void setUp( const GetPot& dataFile );

    void initialize( const Function& d0, const Function& w0 );
    void initialize( const vector_type&, const vector_type& );

    void postProcess();

    void resetPrec() {M_resetPrec = true;}

    //const Dof& dDof() const { return M_FESpace.dof(); }

    //const Mesh& mesh() const { return M_FESpace.mesh(); }
    void reduceSolution( Vector& disp,
                         Vector& vel );

    //Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.getRepeatedEpetra_Map(); }
    EpetraMap const& getMap() const { return M_localMap; }

    const Epetra_Comm& comm() const {return *M_comm;}

    bool isLeader() const
    {
        return comm().MyPID() == 0;
    }

    void evalResidual(bchandler_raw_type & bcFluid, bchandler_raw_type & bcSolid, const vector_type& sol, vector_type& res/*, matrix_type& bigMatrix*/); // used for monolithic

    void evalResidual( const vector_type& sol, vector_type& res); // used for monolithic

    void setDispSolid(const vector_type disp){M_dispSolid = disp;} // used for monolithic

    void updateStuff(); // used for monolithic
    void    updateMatrix(matrix_type & bigMatrixStokes); // used for monolithic
    void    rescaleMatrices(); // used for monolithic
    //    matrix_ptrtype getSolidBlockPtr(){return M_massStiff;}// used for monolithic
    void setBlockPreconditioner(matrix_ptrtype blockPrec){*blockPrec += *M_massStiff;}
    void setFullPreconditioner(matrix_ptrtype fullPrec){*fullPrec += *M_matrFull;}
    matrix_ptrtype getMatrixPtr(){return M_matrFull;}// used for monolithic
private:

    const data_type&               M_data;

    FESpace<Mesh, EpetraMap>&      M_FESpace;

    Epetra_Comm*                   M_comm;
    int                            M_me;
    bool                           M_verbose;

    BCHandler*                     M_BCh;

    EpetraMap                      M_localMap;


    //! Matrix M: mass
    matrix_ptrtype                    M_mass;

    //! Matrix Kl: stiffness linear
    matrix_ptrtype                    M_linearStiff;

    //! Matrix Knl: stiffness non-linear
    matrix_ptrtype                    M_stiff;

    //! Matrix C: mass + linear stiffness
    matrix_ptrtype                    M_massStiff;

    //! Matrix J: jacobian
    matrix_type                    M_jacobian;

    //! Elementary matrices and vectors
    ElemMat                        M_elmatK; // stiffnes
    ElemMat                        M_elmatM; // mass
    ElemMat                        M_elmatC; // mass + stiffness
    ElemVec                        M_elvec;  // Elementary right hand side
    ElemVec                        M_dk_loc; // Local displacement

    //! linearized velocity

    vector_type                    M_disp;
    vector_type                    M_vel;

    //! right  hand  side displacement
    vector_type                    M_rhs;

    //! right  hand  side velocity
    vector_type                    M_rhsW;

    //! right  hand  side
    vector_type                    M_rhsNoBC;

    //! right  hand  side
    vector_type                    M_f;

    //! residual
    vector_type                    M_residual_d;

//    vector_type*                   M_sxx;

    vector_type                    M_sxx;
    vector_type                    M_syy;
    vector_type                    M_szz;

    //! files for lists of iterations and residuals per timestep
    std::ofstream                  M_out_iter;
    std::ofstream                  M_out_res;

    //! level of recursion for Aztec (has a sens with FSI coupling)
    UInt _recur;

    source_type                    M_source;

    //! data for solving tangent problem with aztec
    solver_type                    M_linearSolver;
    prec_type                      M_prec;

    bool                           M_reusePrec;
    int                            M_maxIterForReuse;
    bool                           M_resetPrec;

    int                            M_maxIterSolver;

    int                            M_count;

    UInt                           M_offset;
    bool                           M_monolithic;
    vector_type                    M_dispSolid;
    boost::shared_ptr<FESpace<Mesh, EpetraMap> >      M_uFESpace;
    //    bool                           M_isDiagonalBlockPrec;
    //
    //! methods
    //

    matrix_ptrtype                  M_matrFull;

    void applyBoundaryConditions(matrix_type &matrix,
                                 vector_type &rhs,
                                 bchandler_raw_type& BCh,
                                 UInt         offset=0);


    UInt dim() const { return M_FESpace.dim(); }

};


//
//                                         IMPLEMENTATION
//
template <typename Mesh, typename SolverType>
VenantKirchhofSolver<Mesh, SolverType>::
VenantKirchhofSolver( const data_type&          data,
                      FESpace<Mesh, EpetraMap>& dFESpace,
                      BCHandler&                BCh,
                      Epetra_Comm&              comm,
                      UInt                      offset
                      ) :
    M_data                       ( data ),
    M_FESpace                    ( dFESpace ),
    M_BCh                        ( &BCh ),
    M_comm                       ( &comm ),
    M_me                         ( comm.MyPID() ),
    M_verbose                    ( M_me == 0 ),
    M_linearSolver               ( ),
    M_prec                       ( ),//new prec_raw_type() ),
    M_localMap                   ( M_FESpace.map() ),
    M_mass                       ( new matrix_type(M_localMap) ),
    M_linearStiff                ( new matrix_type(M_localMap) ),
    M_stiff                      ( new matrix_type(M_localMap) ),
    M_massStiff                  ( new matrix_type(M_localMap) ),
    M_jacobian                   ( M_localMap ),
    M_elmatK                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatM                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatC                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elvec                      ( M_FESpace.fe().nbNode, nDimensions ),
    M_dk_loc                     ( M_FESpace.fe().nbNode, nDimensions ),
    M_disp                       ( M_localMap ),
    M_vel                        ( M_localMap ),
    M_rhs                        ( M_localMap ),
    M_rhsW                       ( M_localMap ),
    M_rhsNoBC               ( M_localMap ),
    M_f                          ( M_localMap),
    M_residual_d                 ( M_localMap ),
    M_sxx                        ( M_localMap ),
    M_syy                        ( M_localMap ),
    M_szz                        ( M_localMap ),
    M_out_iter                   ( "out_iter_solid" ),
    M_out_res                    ( "out_res_solid" ),
    M_reusePrec                  ( true ),
    M_maxIterForReuse            ( -1 ),
    M_resetPrec                  ( true ),
    M_maxIterSolver              ( -1 ),
    M_count                      ( 0 ),
    M_offset                     (offset),
    M_dispSolid                ( M_localMap ),
    M_uFESpace                 ( )
{
    //    M_BCh->setOffset(M_offset);
}

template <typename Mesh, typename SolverType>
VenantKirchhofSolver<Mesh, SolverType>::
VenantKirchhofSolver( const data_type& data,
                      FESpace<Mesh, EpetraMap>&   dFESpace,
                      Epetra_Comm&     comm,
                      UInt             /*offset*/
                      ) :
    M_data                       ( data ),
    M_FESpace                    ( dFESpace ),
    M_comm                       ( &comm ),
    M_me                         ( comm.MyPID() ),
    M_verbose                    ( M_me == 0 ),
    M_localMap                   ( M_FESpace.map() ),
    M_mass                       ( new matrix_type(M_localMap) ),
    M_linearStiff                ( new matrix_type(M_localMap) ),
    M_stiff                      ( new matrix_type(M_localMap) ),
    M_massStiff                  ( new matrix_type(M_localMap) ),
    M_jacobian                   ( M_localMap ),
    M_elmatK                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatM                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatC                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elvec                      ( M_FESpace.fe().nbNode, nDimensions ),
    M_dk_loc                     ( M_FESpace.fe().nbNode, nDimensions ),
    M_disp                       ( M_localMap ),
    M_vel                        ( M_localMap ),
    M_rhs                        ( M_localMap ),
    M_rhsW                       ( M_localMap ),
    M_rhsNoBC                    ( M_localMap ),
    M_f                          ( M_localMap),
    M_residual_d                 ( M_localMap ),
    M_sxx                        ( M_localMap ),
    M_syy                        ( M_localMap ),
    M_szz                        ( M_localMap ),
    M_out_iter                   ( "out_iter_solid" ),
    M_out_res                    ( "out_res_solid" ),
    M_linearSolver               ( ),
    M_prec                       ( ),//new prec_raw_type() ),
    M_reusePrec                  ( true ),
    M_maxIterForReuse            ( -1 ),
    M_resetPrec                  ( true ),
    M_maxIterSolver              ( -1 ),
    M_count                      ( 0 ),
    M_offset                     ( 0 ),
    M_dispSolid                ( M_localMap ),
    M_uFESpace                 ( )
{
    //    M_BCh->setOffset(0);
}


template <typename Mesh, typename SolverType>
VenantKirchhofSolver<Mesh, SolverType>::
VenantKirchhofSolver( const data_type& data,
                      FESpace<Mesh, EpetraMap>&   dFESpace,
                      Epetra_Comm&     comm,
                      EpetraMap&      monolithicMap,
                      UInt             offset,
                      boost::shared_ptr<FESpace<Mesh, EpetraMap> >   uFESpace
                      ):
    M_data                       ( data ),
    M_FESpace                    ( dFESpace ),
    M_comm                       ( &comm ),
    M_me                         ( comm.MyPID() ),
    M_verbose                    ( M_me == 0 ),
    M_linearSolver               ( ),
    M_prec                       ( ),//new prec_raw_type() ),
    M_elmatK                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatM                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatC                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elvec                      ( M_FESpace.fe().nbNode, nDimensions ),
    M_dk_loc                     ( M_FESpace.fe().nbNode, nDimensions ),
    M_out_iter                   ( "out_iter_solid" ),
    M_out_res                    ( "out_res_solid" ),
    M_reusePrec                  ( true ),
    M_resetPrec                  ( true ),
    M_maxIterSolver              ( -1 ),
    M_count                      ( 0 ),
    M_massStiff                  ( new matrix_type(monolithicMap) ),
    M_localMap                   ( monolithicMap ),
    M_mass                       ( new matrix_type(M_localMap) ),
    M_linearStiff                ( new matrix_type(M_localMap) ),
    M_stiff                      ( new matrix_type(M_localMap) ),
    M_jacobian                   ( M_localMap ),
    M_disp                       ( M_localMap ),
    M_vel                        ( M_localMap ),
    M_rhs                        ( M_localMap ),
    M_rhsW                       ( M_localMap ),
    M_rhsNoBC                    ( M_localMap ),
    M_f                          ( M_localMap),
    M_residual_d                 ( M_localMap ),
    M_sxx                        ( M_localMap ),
    M_syy                        ( M_localMap ),
    M_szz                        ( M_localMap ),
    M_offset                     ( offset ),
    M_dispSolid                  ( M_localMap ),
    M_uFESpace                   ( (uFESpace) )
{

    }

//


template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
setUp( const GetPot& dataFile )
{
    if (M_verbose)
    {
        std::cout << std::endl;
        std::cout << "S-  Displacement unknowns: " << M_FESpace.dof().numTotalDof() << std::endl;
        std::cout << "S-  Computing mass and linear strain matrices ... ";
    }

    M_monolithic = dataFile("problem/monolithic"   , false );
    M_linearSolver.setDataFromGetPot( dataFile, "solid/solver" );
//    M_linearSolver.setMatrix( M_jacobian );
    if (M_verbose)
        std::cout << "ok." << std::endl;

    M_reusePrec     = dataFile( "solid/prec/reuse", true);
    M_maxIterSolver = dataFile( "solid/solver/max_iter", -1);
    M_maxIterForReuse = dataFile( "solid/solver/max_iter_reuse", M_maxIterSolver*8/10);

    std::string precType = dataFile( "solid/prec/prectype", "Ifpack");

    M_prec.reset( PRECFactory::instance().createObject( precType ) );
    ASSERT(M_prec.get() != 0, "VenantKirchhof : Preconditioner not set");
    //    M_prec               = prec_ptr( PRECFactory::instance().createObject( precType ) );

    M_prec->setDataFromGetPot( dataFile, "solid/prec" );


}


template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::
updateMatrix(matrix_type & bigMatrixStokes)
{
    //    M_massStiff->GlobalAssemble();
    bigMatrixStokes.GlobalAssemble();
    if(M_matrFull)
        M_matrFull.reset(new matrix_type(M_localMap, M_matrFull->getMeanNumEntries()));
    else
        M_matrFull.reset(new matrix_type(M_localMap, M_massStiff->getMeanNumEntries()));
    //    M_matrFull += M_massStiff;
    //    M_matrFull->GlobalAssemble();
    //    M_massStiff.reset(new matrix_type(M_localMap));
    *M_matrFull += *M_massStiff;
    *M_matrFull += bigMatrixStokes;

}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::
rescaleMatrices()
{
    *M_mass *=M_data.timestep();
    *M_linearStiff *= M_data.timestep();
    *M_massStiff *= M_data.timestep();
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::
buildSystem()
{
    UInt totalDof = M_FESpace.dof().numTotalDof();

    if (M_verbose)
        std::cout << "S-  Building the system                       ... ";

    Chrono chrono;
    chrono.start();

    // Number of displacement components
    UInt nc = nDimensions;

    //inverse of dt:
    Real dti2 = 2.0 / ( M_data.timestep() * M_data.timestep() );

    // Elementary computation and matrix assembling
    // Loop on elements
    for ( UInt i = 1; i <= M_FESpace.mesh()->numVolumes(); i++ )
    {

        M_FESpace.fe().updateFirstDerivQuadPt( M_FESpace.mesh()->volumeList( i ) );

        // int marker    = M_FESpace.mesh()->volumeList( i ).marker();

        double lambda = M_data.lambda();
        double mu     = M_data.mu    ();

        M_elmatK.zero();
        M_elmatM.zero();


        // stiffness
        stiff_strain(     mu, M_elmatK, M_FESpace.fe() );
        stiff_div   ( 0.5*lambda, M_elmatK, M_FESpace.fe() );

        M_elmatC.mat() = M_elmatK.mat();

        // mass
        mass( dti2 * M_data.rho(), M_elmatM, M_FESpace.fe(), 0, 0, nDimensions );

        M_elmatC.mat() += M_elmatM.mat();

        // assembling
        for ( UInt ic = 0; ic < nc; ic++ )
        {
            for ( UInt jc = 0; jc < nc; jc++ )
            {
                assembleMatrix( *M_linearStiff, M_elmatK, M_FESpace.fe(), M_FESpace.fe(), M_FESpace.dof(), M_FESpace.dof(),  ic,  jc,  M_offset +ic*totalDof, M_offset + jc*totalDof );

                assembleMatrix( *M_massStiff  , M_elmatC, M_FESpace.fe(), M_FESpace.fe(), M_FESpace.dof(), M_FESpace.dof(),  ic,jc, M_offset +ic*totalDof,M_offset + jc*totalDof );
            }

            //mass
            assembleMatrix( *M_mass, M_elmatM, M_FESpace.fe(), M_FESpace.dof(),ic,ic, M_offset +  ic*totalDof, M_offset +  ic*totalDof);
        }
    }

    M_comm->Barrier();

    M_linearStiff->GlobalAssemble();
    M_massStiff->GlobalAssemble();
    M_mass->GlobalAssemble();


//      M_linearStiff->spy("linearStiff");
//      M_massStiff->spy("massStiff");
//      M_mass->spy("mass");
    chrono.stop();

    if (M_verbose)
        std::cout << " solid done in " << chrono.diff() << " s." << std::endl;


}

//


template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
initialize( const vector_type& disp, const vector_type& vel)
{
    M_disp = disp;
    M_vel  = vel;
}



//
template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
initialize( const Function& d0, const Function& w0 )
{

    M_FESpace.interpolate(d0, M_disp, 0.0);
    M_FESpace.interpolate(w0, M_vel , 0.0);

}


template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
updateSystem( vector_type & rhsFluidCoupling )
{


    if (M_verbose)
    {
        std::cout << "  s-  Updating mass term on right hand side... " << std::flush;
    }

    Chrono chrono;
    chrono.start();

    M_stiff = M_linearStiff;

    // Number of displacement components
    UInt nc = nDimensions;

    double coef;

    coef = (double) M_data.timestep();

    //M_dispSolid already initiaalized
    vector_type _z = M_dispSolid;
    _z            += coef*M_vel;

    M_rhsNoBC = rhsFluidCoupling;
    M_rhsNoBC += *M_mass * _z ;
    M_rhsNoBC -= *M_stiff*M_dispSolid;

    coef = 2.0/M_data.timestep();

    M_rhsW  = coef*M_dispSolid;
    M_rhsW += M_vel;

    if (M_verbose) std::cout << std::endl;

//     std::cout << "rhsWithoutBC norm = " << M_rhsNoBC.NormInf() << std::endl;
//     std::cout << "M_rhsW norm       = " << M_rhsW.NormInf() << std::endl;
//     std::cout << "    _w norm       = " << M_vel.NormInf() << std::endl;

    //
    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
}


template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
updateSystem( )
{


    if (M_verbose)
    {
        std::cout << "  s-  Updating mass term on right hand side... " << std::flush;
    }

    Chrono chrono;
    chrono.start();

    M_stiff = M_linearStiff;

    // Number of displacement components
//    UInt nc = nDimensions;

    double coef;

    coef = (double) M_data.timestep();

    vector_type _z = M_disp;
    _z            += coef*M_vel;

    M_rhsNoBC  = *M_mass*_z;
    M_rhsNoBC -= *M_stiff*M_disp;

    coef = 2.0/M_data.timestep();

    M_rhsW  = coef*M_disp;
    M_rhsW += M_vel;

    if (M_verbose) std::cout << std::endl;

//     std::cout << "rhsWithoutBC norm = " << M_rhsNoBC.NormInf() << std::endl;
//     std::cout << "M_rhsW norm       = " << M_rhsW.NormInf() << std::endl;
//     std::cout << "    _w norm       = " << M_vel.NormInf() << std::endl;

    //
    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
}



template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
iterate( bchandler_raw_type& bch )
{
    Chrono chrono;

    std::cout << "  S-  Solving the system ... \n"  << std::flush;

    std::cout << "  S-  Updating the boundary conditions ... " << std::flush;

    chrono.start();


    matrix_ptrtype matrFull( new matrix_type( M_localMap, M_massStiff->getMeanNumEntries()));
    *matrFull += *M_massStiff;

    M_rhsNoBC.GlobalAssemble();
    M_rhsW.GlobalAssemble();

    vector_type rhs (M_rhsNoBC);

    applyBoundaryConditions( *matrFull, rhs, bch );

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s \n" <<  std::flush;

//    M_comm->Barrier();

    M_linearSolver.setMatrix(*matrFull);

    if ( !M_reusePrec || M_resetPrec || !M_prec->set() )
    {
        chrono.start();

        if (M_verbose)
            std::cout << "  S-  Computing the precond ...                " <<  std::flush;

        M_prec->buildPreconditioner(matrFull);

//    M_disp *= 0.;
        double condest = M_prec->Condest();

        M_linearSolver.setPreconditioner(M_prec);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "         Estimated condition number = " << condest << "\n" <<  std::flush;
        }

        M_resetPrec = false;

    }
    else
    {
        if (M_verbose)
            std::cout << "  S-  Reusing  precond ...                \n" <<  std::flush;
    }

    int numIter = M_linearSolver.solve(M_disp, rhs);

    M_vel  = ( 2.0 / M_data.timestep() ) * M_disp;
    M_vel -= M_rhsW;

    M_residual_d =  *M_massStiff * M_disp;
    M_residual_d -= M_rhsNoBC;
//    reduceSolution(M_disp, M_vel);

    std::cout << "  . " << std::endl;
    if (M_verbose)
    {
        std::cout << "S- system solved in " << chrono.diff()
                  << " s. ( " << numIter << "  iterations. ) \n"
                  << std::flush;
    }


    if (numIter > M_maxIterSolver)
    {
        M_resetPrec = true;
    }



//    postProcess();

}

template <typename Mesh, typename SolverType> // for monolithic
void VenantKirchhofSolver<Mesh, SolverType>::
evalResidual( const vector_type& sol,vector_type& res )
{
    res = *M_matrFull*sol;
    res -= M_rhs;
    //    res.spy("res");
    // Ax-b
}
template <typename Mesh, typename SolverType> // for monolithic
void VenantKirchhofSolver<Mesh, SolverType>::
evalResidual( bchandler_raw_type& bchFluid, bchandler_raw_type& bchSolid, const vector_type& sol,vector_type& res/*, matrix_type& bigMatrixStokes*/ )
{

    //    M_matrFull.reset( new matrix_type(*M_massStiff) );

    vector_type rhs (M_rhsNoBC);

    //    rhs.spy("rhs0");

    applyBoundaryConditions( *M_matrFull, rhs, bchSolid, M_offset );


    //    applyBoundaryConditions( *M_matrFull, rhs, bchFluid, 0 );
    if ( !bchFluid.bdUpdateDone() )
        bchFluid.bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );

    //  vector_type rhsFull(rhs, Repeated/*Unique?*/, Zero); // ignoring non-local entries, Otherwise they are summed up lately
    vector_type rhsFull(rhs, Unique);  // Why there is a ths, rhsFull, and a M_rhs, in this method?

    bcManage( *M_matrFull, rhsFull, *M_uFESpace->mesh(), M_uFESpace->dof(), bchFluid, M_uFESpace->feBd(), 1.,
              M_data.time() );

    std::cout << " rhs=> " << rhs.NormInf() << std::endl;

    //    rhs.spy("rhs2");
    M_rhs=rhsFull;

    res = *M_matrFull*sol;
    res -= M_rhs;
    //    res.spy("res");
    //    M_rhs=res;
}



template <typename Mesh, typename SolverType> // for monolithic
void VenantKirchhofSolver<Mesh, SolverType>::
iterateMonolithic(vector_type& rhs, vector_type& step, matrix_ptrtype prec)
{
    Chrono chrono;

    std::cout << "  S-  Solving the system ... " << std::endl << std::flush;

    std::cout << "  S-  Updating the boundary conditions ... " << std::flush;

    chrono.start();

    // boundary conditions applied in the residual evaluation

    M_matrFull->spy("jacobian");
    prec->spy("blockPreconditioner");
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s \n";

    /*        for (UInt i = 7980; i < 9240; i++ )//lines to kill
        {
            double* tmp;
            int err;
            int entries;
            err = M_matrFull->getEpetraMatrix().ExtractGlobalRowView(i, entries, tmp);
            if(entries == 0)
                {
                    std::cout<<"ERROR in line " << i << " err " << err << std::endl;
                    //                    break;
                }
            //            std::cout << "nonzero entries for row " << i << " ==> " << entries << std::endl;
            }*/

    M_comm->Barrier();

    M_linearSolver.setMatrix(*M_matrFull);

    if  ( !M_reusePrec || M_resetPrec )
    {
        chrono.start();

        if (M_verbose)
            std::cout << "  S-  Computing the precond ...                ";

        M_prec->buildPreconditioner(prec/*M_matrFull*/);

//    M_disp *= 0.;
        double condest = M_prec->Condest();

        M_linearSolver.setPreconditioner(M_prec);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "         Estimated condition number = " << condest << "\n" <<  std::flush;
        }

        M_resetPrec = false;

    }
    else
    {
        if (M_verbose)
            std::cout << "  S-  Reusing  precond ...                \n" <<  std::flush;
    }
    //    M_disp.spy("disp0");
    int numIter = M_linearSolver.solve(step, rhs);

    if (numIter > M_maxIterSolver)
    {
        chrono.start();

        std::cout<<"  s- Iterative solver failed, numiter = " << numIter<<std::endl;

        M_prec->buildPreconditioner(M_matrFull);

        double condest = M_prec->Condest();

        M_linearSolver.setPreconditioner(M_prec);

        chrono.stop();

        numIter = M_linearSolver.solve(step, rhs);

        if (numIter > M_maxIterSolver && M_verbose)
            std::cout << "  s- ERROR: Iterative solver failed again.\n" <<  std::flush;
    }

    M_disp += step;

    std::cout << "  S- system solved. " << std::endl;

    //    M_dispSolid.spy("dispSolid0");
}
template <typename Mesh, typename SolverType> // for monolithic
void VenantKirchhofSolver<Mesh, SolverType>::
updateStuff()
{
    M_vel  = ( 2.0 / M_data.timestep() ) * M_dispSolid;
    M_vel -= M_rhsW;

    M_residual_d =  *M_matrFull/**M_massStiff*/*M_disp;
    M_residual_d -= M_rhsNoBC;
//    reduceSolution(M_disp, M_vel);
//    M_dispSolid.spy("dispSolid");
//    M_disp.spy("disp");


//    postProcess();

}



template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
iterateLin( bchandler_raw_type& bch )
{
    Chrono chrono;

    std::cout << "  S-  Solving the system ... " << std::endl << std::flush;

    std::cout << "  S-  Updating the boundary conditions ... " << std::flush;

    chrono.start();

    matrix_ptrtype matrFull( new matrix_type( M_localMap, M_massStiff->getMeanNumEntries()));
    *matrFull += *M_massStiff;

    vector_type rhs(M_FESpace.map());
    rhs *= 0.;

    applyBoundaryConditions( *matrFull, rhs, bch );

    std::cout << "rhs_dz norm = " << rhs.NormInf() << std::endl;

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s \n";

//    M_comm->Barrier();

    M_linearSolver.setMatrix(*matrFull);

    if  ( !M_reusePrec || M_resetPrec )
    {
        chrono.start();

        if (M_verbose)
            std::cout << "  S-  Computing the precond ...                ";

        M_prec->buildPreconditioner(matrFull);

//    M_disp *= 0.;
        double condest = M_prec->Condest();

        M_linearSolver.setPreconditioner(M_prec);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "         Estimated condition number = " << condest << "\n" <<  std::flush;
        }

        M_resetPrec = false;

    }
    else
    {
        if (M_verbose)
            std::cout << "  S-  Reusing  precond ...                \n" <<  std::flush;
    }

    int numIter = M_linearSolver.solve(M_disp, rhs);

    std::cout << "dz norm     = " << M_disp.NormInf() << std::endl;


    M_residual_d =  *M_massStiff*M_disp;
//    M_residual_d -= M_rhsNoBC;
//    reduceSolution(M_disp, M_vel);

    std::cout << "  S- system solved. " << std::endl;

    if (numIter > M_maxIterSolver)
    {
        M_resetPrec = true;
    }



//    postProcess();

}













template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
iterate(vector_type &_sol)
{

    int status;

    vector_type sol(M_localMap);

    M_vel  = ( 2.0 / M_data.timestep() ) * M_disp;
    M_vel -= M_rhsW;

    std::cout << "sol norm = " << norm(this->sol) << std::endl;

    M_residual_d  = M_massStiff*sol;
    M_residual_d -= M_rhsNoBC;
}


template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
evalResidual( vector_type &res, const vector_type& sol, int /*iter*/)
{
    std::cout << "    s-    Computing residual... ";
    Chrono chrono;
    chrono.start();

    // Matrices initialization
    M_stiff = M_massStiff;


    std::cout << "updating the boundary conditions" << std::flush;
    if ( !M_BCh->bdUpdateDone() )
        M_BCh->bdUpdate( M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );
    std::cout << std::endl;

    bcManageMatrix( M_stiff, *M_FESpace.mesh(), M_FESpace.dof(), *M_BCh, M_FESpace.feBd(), 1.0 );

    M_rhs = M_rhsNoBC;

    bcManageVector( M_rhs, *M_FESpace.mesh(), M_FESpace.dof(), *M_BCh, M_FESpace.feBd(), M_data.time(), 1.0 );

    res  = M_stiff * sol;
//    res -= M_rhs;

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
}



template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
updateJacobian( vector_type & sol, int iter )
{
    std::cout << "  S-  Solid: Updating JACOBIAN in iter " << iter << "  ... ";

    Chrono chrono;
    chrono.start();


    // copy of the linear part
    M_jacobian = M_massStiff;

    /*
    if ( _maxiter > 1 )
    {
        std::cout << " ******** non linear part" << std::endl;

        UInt ig;

        // Number of displacement components
        UInt nc = this->_d.nbcomp();

        // loop on volumes: assembling source term
        for ( UInt i = 1; i <= M_FESpace.mesh()->numVolumes(); ++i )
        {
            M_FESpace.fe().updateFirstDerivQuadPt( M_FESpace.mesh()->volumeList( i ) );

            _elmatK.zero();

            // _dk_loc contains the displacement in the nodes
            for ( UInt j = 0 ; j < ( UInt ) M_FESpace.fe().nbNode ; ++j )
            {
                for ( UInt ic = 0; ic < nc; ++ic )
                {
                    ig = M_FESpace.dof().localToGlobal( i, j + 1 ) - 1 + ic * dim();
                    _dk_loc[ j + ic * M_FESpace.fe().nbNode ] = sol[ ig ];
                }
            }

            // stiffness for non-linear terms
            // 1/2 * \mu * ( [\grad \delta d]^T \grad d^k + [\grad d^k]^T \grad \delta d : \grad v  )
            stiff_dergrad( this->_mu * 0.5, _dk_loc, _elmatK, M_FESpace.fe() );

            // 1/2 * \lambda * ( \tr { [\grad u^k]^T \grad u }, \div v  )
            stiff_derdiv( 0.5 * this->_lambda, _dk_loc, _elmatK, M_FESpace.fe() );

            // assembleing
            for ( UInt ic = 0; ic < nc; ++ic )
                for ( UInt jc = 0; jc < nc; jc++ )
                    assemb_mat( M_jacobian, _elmatK, M_FESpace.fe(), M_FESpace.dof(), ic, jc );
        }
    }
    */
//     if (iter == 1)
//     {
//         M_jacobian.spy("Jacobian");
//     }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
}



template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
evalConstraintTensor()
{

//    double tmp;

    vector_type count(M_localMap);

    M_sxx *= 0.;
    M_syy *= 0.;
    M_szz *= 0.;

   for ( UInt ielem = 1; ielem <= M_FESpace.mesh()->numVolumes(); ielem++ )
    {
        //UInt elem = M_FESpace.mesh()->volumeList( ielem ).id();
        M_FESpace.fe().updateFirstDerivQuadPt( M_FESpace.mesh()->volumeList( ielem ) );

        //int    marker = M_FESpace.mesh()->volumeList( ielem ).marker();
        double s      = 0;
        double volume = M_FESpace.fe().detJac(0);

        for ( int ig = 0; ig < M_FESpace.fe().nbQuadPt; ++ig )
        {
            for ( int k = 0; k < M_FESpace.fe().nbNode; ++k )
            {
                int i    = M_FESpace.fe().patternFirst(k);
                int idof = M_FESpace.dof().localToGlobal(M_FESpace.fe().currentLocalId(), i + 1);

                s+= (2*M_data.mu() + M_data.lambda())*
                    M_FESpace.fe().weightDet( ig )*
                    M_FESpace.fe().phiDer( k, 0 , ig )*
                    M_disp[idof + 0*M_FESpace.dim()];

                s+= M_data.lambda()*
                    M_FESpace.fe().weightDet( ig )*
                    M_FESpace.fe().phiDer( k, 1 , ig )*
                    M_disp[idof + 1*M_FESpace.dim()];

                s+= M_data.lambda()*
                    M_FESpace.fe().weightDet( ig )*
                    M_FESpace.fe().phiDer( k, 2 , ig )*
                    M_disp[idof + 2*M_FESpace.dim()];

                count[idof]++;
            }
        }

        for ( int k = 0; k < M_FESpace.fe().nbNode; ++k )
        {
            int i    = M_FESpace.fe().patternFirst(k);
            int idof = M_FESpace.dof().localToGlobal(M_FESpace.fe().currentLocalId(), i + 1);

            M_sxx[idof] += s/M_FESpace.fe().detJac(0);
        }

        s = 0;

        for ( int ig = 0; ig < M_FESpace.fe().nbQuadPt; ++ig )
        {
            for ( int k = 0; k < M_FESpace.fe().nbNode; ++k )
            {
                int i    = M_FESpace.fe().patternFirst(k);
                int idof = M_FESpace.dof().localToGlobal(M_FESpace.fe().currentLocalId(), i + 1);

                s += M_data.lambda()*
                    M_FESpace.fe().weightDet( ig )*
                    M_FESpace.fe().phiDer( k, 0 , ig )*
                    M_disp[idof + 0*M_FESpace.dim()];

                s += (2*M_data.mu() + M_data.lambda())*
                    M_FESpace.fe().weightDet( ig )*
                    M_FESpace.fe().phiDer( k, 1 , ig )*
                    M_disp[idof + 1*M_FESpace.dim()];

                s += M_data.lambda()*
                    M_FESpace.fe().weightDet( ig )*
                    M_FESpace.fe().phiDer( k, 2 , ig )*
                    M_disp[idof + 2*M_FESpace.dim()];
//                         M_sxx[idof] += s;

//                M_syy[idof] += s/volume;
            }
        }

         for ( int k = 0; k < M_FESpace.fe().nbNode; ++k )
         {
             int i    = M_FESpace.fe().patternFirst(k);
             int idof = M_FESpace.dof().localToGlobal(M_FESpace.fe().currentLocalId(), i + 1);

             M_syy[idof] += s/volume;
         }


         s = 0;

        for ( int ig = 0; ig < M_FESpace.fe().nbQuadPt; ++ig )
        {
            for ( int k = 0; k < M_FESpace.fe().nbNode; ++k )
            {
                int i    = M_FESpace.fe().patternFirst(k);
                int idof = M_FESpace.dof().localToGlobal(M_FESpace.fe().currentLocalId(), i + 1);

                s += M_data.lambda()*
                    M_FESpace.fe().weightDet( ig )*
                    M_FESpace.fe().phiDer( k, 0 , ig )*
                    M_disp[idof + 0*M_FESpace.dim()];

                s += M_data.lambda()*
                    M_FESpace.fe().weightDet( ig )*
                    M_FESpace.fe().phiDer( k, 1 , ig )*
                    M_disp[idof + 1*M_FESpace.dim()];

                s += (2*M_data.mu() + M_data.lambda())*
                    M_FESpace.fe().weightDet( ig )*
                    M_FESpace.fe().phiDer( k, 2 , ig )*
                    M_disp[idof + 2*M_FESpace.dim()];


//                         M_sxx[idof] += s;
            }
        }

        for ( int k = 0; k < M_FESpace.fe().nbNode; ++k )
        {
            int i    = M_FESpace.fe().patternFirst(k);
            int idof = M_FESpace.dof().localToGlobal(M_FESpace.fe().currentLocalId(), i + 1);

            M_szz[idof] += s/M_FESpace.fe().detJac(0);
        }

    }

    for (int ii = 1; ii <= (int)M_FESpace.dim(); ++ii)
    {
        M_sxx[ii] /= count[ii];
        M_syy[ii] /= count[ii];
        M_szz[ii] /= count[ii];
    }


}

//     for ( int k = 0; k < M_FESpace.fe().nbNode; ++k )
//     {
//         int i    = M_FESpace.fe().patternFirst(k);
//         int idof = M_FESpace.dof().localToGlobal(M_FESpace.fe().currentLocalId(), i + 1);

//         M_sxx[idof] += sxx/M_FESpace.fe().detJac(0);
//         M_syy[idof] += syy/M_FESpace.fe().detJac(0);
//         M_szz[idof] += szz/M_FESpace.fe().detJac(0);

//     }

//     for (int ii = 1; ii <= M_FESpace.dim(); ++ii)
//     {
//         M_sxx[ii] /= count[ii];
//         M_syy[ii] /= count[ii];
//         M_szz[ii] /= count[ii];
//     }




template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
//solveJac( const Vector& res, double& linear_rel_tol, Vector &step)
solveJac( vector_type &step, const vector_type& res, double& /*linear_rel_tol*/)
{
    Chrono chrono;

    M_f = res;

    // for BC treatment (done at each time-step)
    Real tgv = 1.0;
    std::cout << "   S-  Applying boundary conditions      ... ";
    chrono.start();

    // BC manage for the velocity
    if ( !M_BCh->bdUpdateDone() )
        M_BCh->bdUpdate( M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );

    bcManageMatrix( M_jacobian, *M_FESpace.mesh(), M_FESpace.dof(), *M_BCh, M_FESpace.feBd(), tgv );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    M_jacobian.spy("jacobian");
//    M_linearSolver.setRecursionLevel( _recur );

    std::cout << "   S-  Solving system                    ... "<< std::flush;
    chrono.start();
    M_linearSolver.solve( step , M_f);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    //--options[AZ_recursion_level];

//    AZ_matrix_destroy( &J );
//    AZ_precond_destroy( &precM_jacobian );

    M_residual_d = M_massStiff*step;// - M_rhsNoBC;
}


template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
//solveJac( const Vector& res, double& linear_rel_tol, Vector &step)
solveJac( vector_type&       step,
          const vector_type& res,
          double&            /*linear_rel_tol*/,
          bchandler_type&    BCh)
{
    Chrono chrono;

    M_f = res;

    // for BC treatment (done at each time-step)
    Real tgv = 1.0;

    std::cout << "   S-  Applying boundary conditions      ... ";
    chrono.start();

    // BC manage for the velocity
    if ( BCh->bdUpdateDone() )
        BCh->bdUpdate( M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );

    bcManageMatrix( M_jacobian, *M_FESpace.mesh(), M_FESpace.dof(), *BCh, M_FESpace.feBd(), tgv );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

//    M_linearSolver.setRecursionLevel( _recur );

    std::cout << "   S-  Solving system                    ... "<< std::flush;
    chrono.start();
//    M_linearSolver.solve( step , _f);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    //--options[AZ_recursion_level];

//    AZ_matrix_destroy( &J );
//    AZ_precond_destroy( &prec_J );

    M_residual_d = M_massStiff*step;// - M_rhsNoBC;
}


template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
solveJacobian( Real /*time*/ )
{

    std::cout << "  S-  LINEARIZED SOLID SYSTEM" << std::flush << std::endl;
    Chrono chrono;

    //if (BCd == 0) BCd.reset(&BCh_solid());

//    _f = ZeroVector( _f.size() );
    M_jacobian = M_massStiff;

    // for BC treatment (done at each time-step)
    double tgv = 1.0;

    std::cout << "TGV = " << tgv << std::flush << std::endl;
    std::cout << "  S-  Applying boundary conditions        ... ";
    chrono.start();

    // BC manage for the velocity
    if ( !M_BCh->bdUpdateDone() )
        M_BCh->bdUpdate( M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );

    bcManageVector(M_f,
                   *M_FESpace.mesh(),
                   M_FESpace.dof(),
                   M_BCh,
                   M_FESpace.feBd(),
                   1., 1.);

    bcManageMatrix( M_jacobian, *M_FESpace.mesh(), M_FESpace.dof(), M_BCh, M_FESpace.feBd(), tgv );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

//    M_linearSolver.setRecursionLevel( _recur );

    std::cout << "  S-  Solving system                      ... " << std::flush;
    chrono.start();

    M_linearSolver.solve( M_disp , M_f );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;


    this->_w = ( 2.0 / M_data.timestep() ) * M_disp - M_rhsW;

    M_residual_d = M_massStiff*M_disp;
}

template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
solveJacobian( const Real /*time*/ , bchandler_type& BCd)
{
    std::cout << "  S-  LINEARIZED SOLID SYSTEM" << std::flush << std::endl;
    Chrono chrono;

    //if (BCd == 0) BCd.reset(&BCh_solid());

//    _f = ZeroVector( _f.size() );
    M_jacobian = M_massStiff;

    // for BC treatment (done at each time-step)
    Real tgv = 1.0;

    std::cout << "  S-  Applying boundary conditions        ... ";
    chrono.start();

    // BC manage for the velocity
    if ( !(*BCd).bdUpdateDone() )
        (*BCd).bdUpdate( M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );

    bcManageVector(M_f,
                   *M_FESpace.mesh(),
                   M_FESpace.dof(),
                   *BCd,
                   M_FESpace.feBd(),
                   1., 1.);

    bcManageMatrix( M_jacobian, *M_FESpace.mesh(), M_FESpace.dof(), *BCd, M_FESpace.feBd(), tgv );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

//    M_linearSolver.setRecursionLevel( _recur );

    std::cout << "  S-  Solving system                      ... "<< std::flush;

    chrono.start();
//    M_linearSolver.solve( M_disp , _f );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    this->_w = ( 2.0 / M_data.timestep() ) * M_disp - M_rhsW;

    M_residual_d = M_massStiff*M_disp;
}

template<typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
applyBoundaryConditions(matrix_type&        matrix,
                        vector_type&        rhs,
                        bchandler_raw_type& BCh,
                        UInt                offset)
{

    // BC manage for the velocity
    if(offset)
        BCh.setOffset(offset);
    if ( !BCh.bdUpdateDone() )
        BCh.bdUpdate( *M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );

    // vector_type rhsFull(rhs, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately
    vector_type rhsFull(rhs, Unique);  // bcManages now manages the also repeated parts

    bcManage( matrix, rhsFull, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 1.,
              M_data.time() );

    // matrix should be GlobalAssembled by  bcManage

    rhs = rhsFull;

} // applyBoundaryCondition


template<typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
reduceSolution( Vector& disp,
                Vector& vel)
{
    vector_type displacement(M_disp, 0);
    vector_type velocity    (M_vel , 0);

    if (M_verbose)
    {
            for ( UInt iDof = 0; iDof < nDimensions*dim(); ++iDof )
            {
                disp[ iDof ] = displacement[ iDof + 1 ];
                vel [ iDof ] = velocity    [ iDof + 1 ];
            }
    }
}

template<typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::
VenantKirchhofSolver<Mesh, SolverType>::postProcess()
{
    std::ostringstream index, indexMe;
    std::string name, me, namedef;

    PhysVectUnknown<Vector> disp(dim());
    PhysVectUnknown<Vector> vel (dim());

    reduceSolution(disp, vel);
    evalConstraintTensor();

    M_count++;

//    std::cout << "factor " << M_data.factor() << std::endl;
     if (M_verbose)
     {
        if ( fmod( float( M_count ), float( M_data.verbose() ) ) == 0.0 )
        {
            if (M_verbose) std::cout << "  S-  Post-processing \n";
            index << ( M_count / M_data.verbose() );

            switch ( index.str().size() )
            {
                case 1:
                    name = "00" + index.str();
                    break;
                case 2:
                    name = "0" + index.str();
                    break;
                case 3:
                    name = index.str();
                    break;
            }

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

            vector_type disp(M_disp, Repeated);
            vector_type vel (M_vel,  Repeated);


            namedef = "defor." + me + "." + name + ".mesh";

            meditSolutionWriter( "dep_x." + me + "." + name + ".bb", *M_FESpace.mesh(), disp, M_FESpace.dof().numTotalDof()*0);
            meditSolutionWriter( "dep_y." + me + "." + name + ".bb", *M_FESpace.mesh(), disp, M_FESpace.dof().numTotalDof()*1);
            meditSolutionWriter( "dep_z." + me + "." + name + ".bb", *M_FESpace.mesh(), disp, M_FESpace.dof().numTotalDof()*2);

//             meditSolutionWriter( "res_x." + me + "." + name + ".bb", *M_FESpace.mesh(), M_sxx, M_FESpace.dof().numTotalDof()*0);
//             meditSolutionWriter( "res_y." + me + "." + name + ".bb", *M_FESpace.mesh(), M_syy, M_FESpace.dof().numTotalDof()*1);
//             meditSolutionWriter( "res_z." + me + "." + name + ".bb", *M_FESpace.mesh(), M_szz, M_FESpace.dof().numTotalDof()*2);

            meditSolutionWriter( "resd_x." + me + "." + name + ".bb", *M_FESpace.mesh(), M_residual_d, M_FESpace.dof().numTotalDof()*0);
            meditSolutionWriter( "resd_y." + me + "." + name + ".bb", *M_FESpace.mesh(), M_residual_d, M_FESpace.dof().numTotalDof()*1);
            meditSolutionWriter( "resd_z." + me + "." + name + ".bb", *M_FESpace.mesh(), M_residual_d, M_FESpace.dof().numTotalDof()*2);

//             wr_medit_ascii_scalar( "dep_x." + name + "." + me + ".bb", disp.getEpetraVector().Values() + 0*dim(), M_FESpace.mesh()->numVertices() );
//             wr_medit_ascii_scalar( "dep_y." + name + "." + me + ".bb", disp.getEpetraVector().Values() + 1*dim(), M_FESpace.mesh()->numVertices() );
//             wr_medit_ascii_scalar( "dep_z." + name + "." + me + ".bb", disp.getEpetraVector().Values() + 2*dim(), M_FESpace.mesh()->numVertices() );

            wr_medit_ascii2( namedef, *M_FESpace.mesh(), disp, M_data.factor() );
            // wr_medit_ascii_vector("veloc."+name+".bb",_u.giveVec(),M_FESpace.mesh()->numVertices(),_dim_u);

            system( ( "ln -sf " + namedef + " dep_x." + me + "." + name + ".mesh" ).data() );
            system( ( "ln -sf " + namedef + " dep_y." + me + "." + name + ".mesh" ).data() );
            system( ( "ln -sf " + namedef + " dep_z." + me + "." + name + ".mesh" ).data() );

//             system( ( "ln -sf " + namedef + " res_x." + me + "." + name + ".mesh" ).data() );
//             system( ( "ln -sf " + namedef + " res_y." + me + "." + name + ".mesh" ).data() );
//             system( ( "ln -sf " + namedef + " res_z." + me + "." + name + ".mesh" ).data() );

            system( ( "ln -sf " + namedef + " resd_x." + me + "." + name + ".mesh" ).data() );
            system( ( "ln -sf " + namedef + " resd_y." + me + "." + name + ".mesh" ).data() );
            system( ( "ln -sf " + namedef + " resd_z." + me + "." + name + ".mesh" ).data() );

            // system(("ln -s "+M_FESpace.mesh()_file+" veloc."+name+".mesh").data());

//             wr_medit_ascii_scalar( "veld_x." + name + "." + me + ".bb", vel.getEpetraVector().Values() + 0*dim(), M_FESpace.mesh()->numVertices() );
//             wr_medit_ascii_scalar( "veld_y." + name + "." + me + ".bb", vel.getEpetraVector().Values() + 1*dim(), M_FESpace.mesh()->numVertices() );
//             wr_medit_ascii_scalar( "veld_z." + name + "." + me + ".bb", vel.getEpetraVector().Values() + 2*dim(), M_FESpace.mesh()->numVertices() );

            // wr_medit_ascii_vector("veloc."+name+".bb",_u.giveVec(),M_FESpace.mesh()->numVertices(),_dim_u);

//             system( ( "ln -sf " + namedef + " veld_x." + name + "." + me + ".mesh" ).data() );
//             system( ( "ln -sf " + namedef + " veld_y." + name + "." + me + ".mesh" ).data() );
//             system( ( "ln -sf " + namedef + " veld_z." + name + "." + me + ".mesh" ).data() );

            // system(("ln -s "+M_FESpace.mesh()_file+" veloc."+name+".mesh").data());

        }
     }
}

}
#endif
