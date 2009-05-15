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
  \file NavierStokesAleSolverPC.h
  \author M.A. Fernandez and A. Gauthier
  \date 05/2003
  \version 1.0

  \brief This file contains a NavierStokes ALE solver  class which implements  a semi-implicit scheme with
         an exact factorization. Preconditioning of the  Schur Complement is done by an algebraic
         Chorin-Temam pressure-corrected preconditioner

*/

#ifndef _NAVIERSTOKESALESOLVERPC_HH
#define _NAVIERSTOKESALESOLVERPC_HH

#include <life/lifesolver/NavierStokesAleHandler.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifefem/assembGeneric.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifealg/algebraic_facto.hpp>
#include <life/lifecore/chrono.hpp>
//#include <life/lifealg/dataAztec.hpp>

namespace LifeV
{
/*!
  \class NavierStokesALESolverPC

   This class implements an NavierStokes ALE solver via exact factorization. Preconditioning of the
   Schur Complement is done by an algebraic Chorin-Temam pressure-corrected preconditioner

*/
template <typename Mesh>
class NavierStokesAleSolverPC:
            public NavierStokesAleHandler<Mesh>
{

public:

    typedef typename NavierStokesHandler<Mesh>::Function Function;
    typedef typename NavierStokesHandler<Mesh>::source_type source_type;

    //! Constructors
    /*!
      \param data_file GetPot data file
      \param refFE_u reference FE for the velocity
      \param refFE_p reference FE for the pressure
      \param Qr_u volumic quadrature rule for the velocity
      \param bdQr_u surface quadrature rule for the velocity
      \param Qr_p volumic quadrature rule for the pressure
      \param bdQr_p surface quadrature rule for the pressure
      \param BCh_u boundary conditions for the velocity
      \param BCh_mesh boundary conditions for the motion harmonic extension
    */
    NavierStokesAleSolverPC( const GetPot& data_file,
                             const RefFE& refFE_u,
                             const RefFE& refFE_p,
                             const QuadRule& Qr_u,
                             const QuadRule& bdQr_u,
                             const QuadRule& Qr_p,
                             const QuadRule& bdQr_p,
                             BCHandler& BCh_u,
                             BCHandler& BCh_mesh );

    NavierStokesAleSolverPC( const GetPot& data_file,
                             const DataNavierStokes<Mesh>& dataNavierStokes,
                             const RefFE& refFE_u,
                             const RefFE& refFE_p,
                             const QuadRule& Qr_u,
                             const QuadRule& bdQr_u,
                             const QuadRule& Qr_p,
                             const QuadRule& bdQr_p,
                             BCHandler& BCh_u,
                             BCHandler& BCh_mesh );

    NavierStokesAleSolverPC( const GetPot& data_file,
                             const RefFE& refFE_u,
                             const RefFE& refFE_p,
                             const QuadRule& Qr_u,
                             const QuadRule& bdQr_u,
                             const QuadRule& Qr_p,
                             const QuadRule& bdQr_p);

//     NavierStokesAleSolverPC( const GetPot& data_file,
//                              const RefFE& refFE_u,
//                              const RefFE& refFE_p,
//                              const QuadRule& Qr_u,
//                              const QuadRule& bdQr_u,
//                              const QuadRule& Qr_p,
//                              const QuadRule& bdQr_p)

    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    void timeAdvance( source_type const& source, const Real& time );

    //! Update convective term, bc treatment and solve the linearized ns system
    void iterate( const Real& time );

//    void iterateTransp( const Real& time );

    void iterateLin( const Real& time, BCHandler& BCh_du );
    void solveJacobian(  const Real& time, BCHandler& BCh_du );

    //LIFEV_DEPREPCATED BCHandler & BC_fluid() {return BCh_HarmonicExtension();}

    Vector& residual();
    Vector  getDeltaLambda() {return this->dt()*_du;}

    void setFullEssential(bool _full)
        {
            _factor_data.setFullEssential(_full);
            _factor_data_jacobian.setFullEssential(_full);
        }
//    void setBC(BCHandler &fluidBC, BCHandler &HamonicExtensionBC);
private:

    //! constructor setup
    void     setUp();

    //! simulation time
    double   M_time;

    //! Block pattern of M_u
    MSRPatt _pattM_u_block;

    //! Pattern for M
    MixedPattern<3, 3, MSRPatt> _pattM_u;

    //! Block pattern of C
    MSRPatt _pattC;

    //! Block pattern of D: Vdiv operator
    CSRPatt _pattD_block;

    //! Pattern for D
    MixedPattern<1, 3, CSRPatt> _pattD;

    //! Block  pattern of trD: transpose Vdiv operator trVdiv
    CSRPatt _pattDtr_block;

    //! Pattern  of trD
    MixedPattern<3, 1, CSRPatt> _pattDtr;

    //! Matrix D: Vdiv operator
    MixedMatr<1, 3, CSRPatt, double> _D;

    //! Matrix trD: transpose Vdiv operator trVdiv
    MixedMatr<3, 1, CSRPatt, double> _trD;

    //! Matrix trD: transpose Vdiv operator trVdiv
    MixedMatr<3, 1, CSRPatt, double> _trDAux;

    //! Matrix HinvDtr:  H^{-1}D^T
    MixedMatr<3, 1, CSRPatt, double> _HinvDtr;

    //! Matrix M_u: Vmass
    MixedMatr<3, 3, MSRPatt, double> _M_u;

    //! Matrix HinvC: H^{-1}C
    MSRMatr<double> _HinvC;

    //! Matrix CStokes
    MSRMatr<double> _CStokes;

    //! Matrix C: 1/dt*Vmass + mu*Vstiff operator + Convective_term
    MSRMatr<double> _C;
    MSRMatr<double> _CAux;


    //! H diag matrix: H= diag( _M_u )/sum( diag( _M_u ) ) where _M_u = mass * rho / dt
    std::vector<double> _H;

    //! Elementary matrices and vectors
    ElemMat _elmatC; //velocity stiffnes
    ElemMat _elmatM_u; //velocity mass
    ElemMat _elmatDtr; // vel_i * pres_j
    ElemVec _elvec; // Elementary right hand side
    ElemVec _elvec_du; // Elementary right hand side for the linearized velocity
    ElemVec _elvec_dp; // Elementary right hand side for the linearized pressure
    ElemVec _w_loc; // Elementary mesh velocity
    ElemVec _uk_loc; // Elementary velocity
    ElemVec _pk_loc; // Elementary pressure
    ElemVec _convect; // Elementary convection velocity
    ElemVec _d_loc; // Elementary displacement for right hand side
    ElemVec _dw_loc; // Elementary mesh velocity for right hand side

    //! The velocity
    PhysVectUnknown<Vector> _un;

    //! The linearized pressure
    ScalUnknown<Vector> _dp;

    //! The linearized velocity
    PhysVectUnknown<Vector> _du;

    //! Right  hand  side for the velocity
    PhysVectUnknown<Vector> _f_u;

    //! Right  hand  side for the velocity
    PhysVectUnknown<Vector> _f_duWithOutBC;

    //! Right  hand  side for pressure
    ScalUnknown<Vector> _f_p;

    //! Right  hand  side for the velocity
    PhysVectUnknown<Vector> _f_uWithOutBC;

    //! The residual of the momentum equations
    PhysVectUnknown<Vector> _residual_u;


    //!  This vector contains the product C^{-1}*trD*P where P is the pressure, solution
    //!  of the system (ii).
    PhysVectUnknown<Vector> _invCtrDP;

    DataAztec _dataAztec_i;
    DataAztec _dataAztec_ii;
    DataAztec _dataAztec_s;

    //! DataFactorisation: data passed to matrix-vector product are stored in the class
    DataFactorisation < MSRMatr<double>, MixedMatr<1, 3, CSRPatt, double>,
    MixedMatr<3, 1, CSRPatt, double>, std::vector<double>, MSRMatr<double>, Vector > _factor_data;

    //! DataFactorisation: data passed to matrix-vector product are stored in the class
    DataFactorisation < MSRMatr<double>, MixedMatr<1, 3, CSRPatt, double>,
    MixedMatr<3, 1, CSRPatt, double>, std::vector<double>, MSRMatr<double>, Vector > _factor_data_jacobian;

};


//
//                                         IMPLEMENTATION
//
template <typename Mesh>
NavierStokesAleSolverPC<Mesh>::
NavierStokesAleSolverPC( const GetPot& data_file, const RefFE& refFE_u, const RefFE& refFE_p, const QuadRule& Qr_u,
                         const QuadRule& bdQr_u, const QuadRule& Qr_p, const QuadRule& bdQr_p, BCHandler& BCh_u,
                         BCHandler& BCh_mesh ) :
        NavierStokesAleHandler<Mesh>( data_file,
                                      refFE_u,
                                      refFE_p,
                                      Qr_u,
                                      bdQr_u,
                                      Qr_p,
                                      bdQr_p,
                                      BCh_u,
                                      BCh_mesh ),
        _pattM_u_block( this->_dof_u ),
        _pattM_u( _pattM_u_block, "diag" ),
        _pattC( this->_dof_u, 3 ),
        _pattD_block( this->_dof_p, this->_dof_u ),
        _pattD( _pattD_block ),
        _pattDtr_block( this->_dof_u, this->_dof_p ),
        _pattDtr( _pattDtr_block ),
        _D( _pattD ),
        _trD( _pattDtr ),
        _trDAux( _pattDtr ),
        _HinvDtr( _pattDtr ),
        _M_u( _pattM_u ),
        _HinvC( _pattC ),
        _CStokes( _pattC ),
        _C( _pattC ),
        _CAux( _pattC ),
        _H( _pattC.nRows() ),
        _elmatC( this->_fe_u.nbNode, nDimensions, nDimensions ),
        _elmatM_u( this->_fe_u.nbNode, nDimensions, nDimensions ),
        _elmatDtr( this->_fe_u.nbNode, nDimensions, 0, this->_fe_p.nbNode, 0, 1 ),
        _elvec( this->_fe_u.nbNode, nDimensions ),
        _elvec_du( this->_fe_u.nbNode, nDimensions ),
        _elvec_dp( this->_fe_p.nbNode, 1 ),
        _w_loc( this->_fe_u.nbNode, nDimensions ),
        _uk_loc( this->_fe_u.nbNode, nDimensions ),
        _pk_loc( this->_fe_p.nbNode, 1 ),
        _convect( this->_fe_u.nbNode, nDimensions ),
        _d_loc( this->_fe_u.nbNode, nDimensions ),
        _dw_loc( this->_fe_u.nbNode, nDimensions ),
        _un( this->_dim_u ),
        _dp( this->_dim_p ),
        _du( this->_dim_u ),
        _f_u( this->_dim_u ),
        _f_duWithOutBC( this->_dim_u ),
        _f_p( this->_dim_p ),
        _f_uWithOutBC( this->_dim_u ),
        _residual_u( this->_dim_u ),
        _invCtrDP( this->_dim_u ),
        _dataAztec_i( data_file, "fluid/aztec_i" ),
        _dataAztec_ii( data_file, "fluid/aztec_ii" ),
        _dataAztec_s( data_file, "fluid/aztec_s" ),
        _factor_data( _C, _D, _trD, _H, _HinvC, _HinvDtr, _invCtrDP, _dataAztec_i, _dataAztec_s, this->bcHandler().hasOnlyEssential(), 1 ),
        _factor_data_jacobian( _C, _D, _trD, _H, _HinvC, _HinvDtr, _invCtrDP, _dataAztec_i, _dataAztec_s, this->bcHandler().hasOnlyEssential(), 2 )
{
    setUp();
}

template <typename Mesh>
NavierStokesAleSolverPC<Mesh>::
NavierStokesAleSolverPC( const GetPot&          data_file,
                         const DataNavierStokes<Mesh>& dataNavierStokes,
                         const RefFE&           refFE_u,
                         const RefFE&           refFE_p,
                         const QuadRule&        Qr_u,
                         const QuadRule&        bdQr_u,
                         const QuadRule&        Qr_p,
                         const QuadRule&        bdQr_p,
                         BCHandler&             BCh_u,
                         BCHandler&             BCh_mesh ):
        NavierStokesAleHandler<Mesh>( data_file,
                                      dataNavierStokes,
                                      refFE_u,
                                      refFE_p,
                                      Qr_u,
                                      bdQr_u,
                                      Qr_p,
                                      bdQr_p,
                                      BCh_u,
                                      BCh_mesh ),
        _pattM_u_block( this->_dof_u ),
        _pattM_u( _pattM_u_block, "diag" ),
        _pattC( this->_dof_u, 3 ),
        _pattD_block( this->_dof_p, this->_dof_u ),
        _pattD( _pattD_block ),
        _pattDtr_block( this->_dof_u, this->_dof_p ),
        _pattDtr( _pattDtr_block ),
        _D( _pattD ),
        _trD( _pattDtr ),
        _trDAux( _pattDtr ),
        _HinvDtr( _pattDtr ),
        _M_u( _pattM_u ),
        _HinvC( _pattC ),
        _CStokes( _pattC ),
        _C( _pattC ),
        _CAux( _pattC ),
        _H( _pattC.nRows() ),
        _elmatC( this->_fe_u.nbNode, nDimensions, nDimensions ),
        _elmatM_u( this->_fe_u.nbNode, nDimensions, nDimensions ),
        _elmatDtr( this->_fe_u.nbNode, nDimensions, 0, this->_fe_p.nbNode, 0, 1 ),
        _elvec( this->_fe_u.nbNode, nDimensions ),
        _elvec_du( this->_fe_u.nbNode, nDimensions ),
        _elvec_dp( this->_fe_p.nbNode, 1 ),
        _w_loc( this->_fe_u.nbNode, nDimensions ),
        _uk_loc( this->_fe_u.nbNode, nDimensions ),
        _pk_loc( this->_fe_p.nbNode, 1 ),
        _convect( this->_fe_u.nbNode, nDimensions ),
        _d_loc( this->_fe_u.nbNode, nDimensions ),
        _dw_loc( this->_fe_u.nbNode, nDimensions ),
        _un( this->_dim_u ),
        _dp( this->_dim_p ),
        _du( this->_dim_u ),
        _f_u( this->_dim_u ),
        _f_duWithOutBC( this->_dim_u ),
        _f_p( this->_dim_p ),
        _f_uWithOutBC( this->_dim_u ),
        _residual_u( this->_dim_u ),
        _invCtrDP( this->_dim_u ),
        _dataAztec_i( data_file, "fluid/aztec_i" ),
        _dataAztec_ii( data_file, "fluid/aztec_ii" ),
        _dataAztec_s( data_file, "fluid/aztec_s" ),
        _factor_data( _C, _D, _trD, _H, _HinvC, _HinvDtr, _invCtrDP, _dataAztec_i, _dataAztec_s, this->bcHandler().hasOnlyEssential(), 1 ),
        _factor_data_jacobian( _C, _D, _trD, _H, _HinvC, _HinvDtr, _invCtrDP, _dataAztec_i, _dataAztec_s, this->bcHandler().hasOnlyEssential(), 2 )
{
    std::cout << "New constructor NavierStokesALESolverPC" << std::endl;
    setUp();
}

template <typename Mesh>
NavierStokesAleSolverPC<Mesh>::
NavierStokesAleSolverPC( const GetPot& data_file,
                         const RefFE& refFE_u,
                         const RefFE& refFE_p,
                         const QuadRule& Qr_u,
                         const QuadRule& bdQr_u,
                         const QuadRule& Qr_p,
                         const QuadRule& bdQr_p ) :
    NavierStokesAleHandler<Mesh>( data_file,
                                  refFE_u,
                                  refFE_p,
                                  Qr_u,
                                  bdQr_u,
                                  Qr_p,
                                  bdQr_p),
    _pattM_u_block       ( this->_dof_u ),
    _pattM_u             ( _pattM_u_block, "diag" ),
    _pattC               ( this->_dof_u, 3 ),
    _pattD_block         ( this->_dof_p, this->_dof_u ),
    _pattD               ( _pattD_block ),
    _pattDtr_block       ( this->_dof_u, this->_dof_p ),
    _pattDtr             ( _pattDtr_block ),
    _D                   ( _pattD ),
    _trD                 ( _pattDtr ),
    _trDAux              ( _pattDtr ),
    _HinvDtr             ( _pattDtr ),
    _M_u                 ( _pattM_u ),
    _HinvC               ( _pattC ),
    _CStokes             ( _pattC ),
    _C                   ( _pattC ),
    _CAux                ( _pattC ),
    _H                   ( _pattC.nRows() ),
    _elmatC              ( this->_fe_u.nbNode, nDimensions, nDimensions ),
    _elmatM_u            ( this->_fe_u.nbNode, nDimensions, nDimensions ),
    _elmatDtr            ( this->_fe_u.nbNode, nDimensions, 0, this->_fe_p.nbNode, 0, 1 ),
    _elvec               ( this->_fe_u.nbNode, nDimensions ),
    _elvec_du            ( this->_fe_u.nbNode, nDimensions ),
    _elvec_dp            ( this->_fe_p.nbNode, 1 ),
    _w_loc               ( this->_fe_u.nbNode, nDimensions ),
    _uk_loc              ( this->_fe_u.nbNode, nDimensions ),
    _pk_loc              ( this->_fe_p.nbNode, 1 ),
    _convect             ( this->_fe_u.nbNode, nDimensions ),
    _d_loc               ( this->_fe_u.nbNode, nDimensions ),
    _dw_loc              ( this->_fe_u.nbNode, nDimensions ),
    _dp                  ( this->_dim_p ),
    _un                  ( this->_dim_u ),
    _du                  ( this->_dim_u ),
    _f_u                 ( this->_dim_u ),
    _f_duWithOutBC       ( this->_dim_u ),
    _f_p                 ( this->_dim_p ),
    _f_uWithOutBC        ( this->_dim_u ),
    _residual_u          ( this->_dim_u ),
    _invCtrDP            ( this->_dim_u ),
    _dataAztec_i         ( data_file, "fluid/aztec_i" ),
    _dataAztec_ii        ( data_file, "fluid/aztec_ii" ),
    _dataAztec_s         ( data_file, "fluid/aztec_s" ),
    _factor_data         ( _C, _D, _trD, _H, _HinvC, _HinvDtr, _invCtrDP, _dataAztec_i, _dataAztec_s, true, 1 ),
    _factor_data_jacobian( _C, _D, _trD, _H, _HinvC, _HinvDtr, _invCtrDP, _dataAztec_i, _dataAztec_s, true, 2 )
{
    setUp();
}



// members



template <typename Mesh>
void NavierStokesAleSolverPC<Mesh>::
setUp()
{
    std::cout << std::endl;
    std::cout << "F-  Pressure unknowns: " << this->_dim_p << std::endl;
    std::cout << "F-  Velocity unknowns: " << this->_dim_u << std::endl << std::endl;
    std::cout << "F-  Computing mass matrix... ";

    Chrono chrono;
    chrono.start();

    // Number of velocity components
    UInt nc_u = this->_u.nbcomp();

    // Initializing mass matrix
    _M_u.zeros();

    Real dti = 1.0 / this->dt();

    // loop on volumes: assembling mass term
    for ( UInt i = 1; i <= this->mesh().numVolumes(); ++i )
    {
        this->_fe_u.updateFirstDerivQuadPt( this->mesh().volumeList( i ) );
        _elmatM_u.zero();
        mass( this->density() * dti, _elmatM_u, this->_fe_u, 0, 0, nc_u );
        for ( UInt ic = 0; ic < nc_u; ++ic )
        {
            assemb_mat( _M_u, _elmatM_u, this->_fe_u, this->_dof_u, ic, ic );
        }
    }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
}


template <typename Mesh>
void NavierStokesAleSolverPC<Mesh>::
timeAdvance( source_type const& source, const Real& time )
{

    std::cout << std::endl;
    std::cout << "F== FLUID: Now we are at time " << time << " s." << std::endl;

    M_time = time;
    // Number of velocity components
    UInt nc_u = this->_u.nbcomp();

    std::cout << "  F-  Updating mass term on right hand side... ";

    Chrono chrono;
    chrono.start();

    // Right hand side for the velocity at time
    _f_uWithOutBC = ZeroVector( _f_uWithOutBC.size() );

    // loop on volumes: assembling source term
    for ( UInt i = 1; i <= this->mesh().numVolumes(); ++i )
    {
        _elvec.zero();
        this->_fe_u.updateFirstDerivQuadPt( this->mesh().volumeList( i ) );
        for ( UInt ic = 0; ic < nc_u; ++ic )
        {
            compute_vec( source, _elvec, this->_fe_u, time, ic ); // compute local vector
            assemb_vec( _f_uWithOutBC, _elvec, this->_fe_u, this->_dof_u, ic ); // assemble local vector into global one
        }
    }
    _f_uWithOutBC += _M_u * this->_u;

    // Save last mesh displacement and fluid velocity
    this->_dispOld = this->harmonicExtension().getDisplacement();
    _un = this->_u;

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
}


template <typename Mesh>
void NavierStokesAleSolverPC<Mesh>::
iterate( const Real& time )
{

    Chrono chrono;

    // Number of velocity components
    UInt nc_u = this->_u.nbcomp();

    std::cout << "  F-  Updating matrices... " << std::flush;

    chrono.start();

    //initialize matrices
    _D.zeros();
    _trD.zeros();
    _M_u.zeros();
    _C.zeros();

    Real dti = 1.0 / this->dt();

    // Loop on elements
    for ( UInt i = 1; i <= this->mesh().numVolumes(); i++ )
    {

        this->_fe_p.update( this->mesh().volumeList( i ) ); // just to provide the id number in the
        // assem_mat_mixed
        this->_fe_u.updateFirstDerivQuadPt( this->mesh().volumeList( i ) );

        // initialization of elementary matrices
        _elmatM_u.zero();
        _elmatDtr.zero();
        _elmatC.zero();

        // mass
        mass( this->density() * dti, _elmatM_u, this->_fe_u, 0, 0, nc_u );

        // stiffness strain
        stiff_strain( 2.0 * this->viscosity(), _elmatC, this->_fe_u );
        _elmatC.mat() += _elmatM_u.mat();

        // Non linear term, Semi-implicit approach
        // u_loc contains the velocity values in the nodes

        for ( UInt k = 0 ; k < ( UInt ) this->_fe_u.nbNode ; k++ )
        {
            UInt iloc = this->_fe_u.patternFirst( k );
            for ( UInt ic = 0; ic < nc_u; ++ic )
            {
                UInt ig = this->_dof_u.localToGlobal( i, iloc + 1 ) - 1 + ic * this->_dim_u;
                _elvec[ iloc + ic * this->_fe_u.nbNode ] = this->density() * ( _un( ig ) - this->_wInterp( ig ) );
                _w_loc[ iloc + ic * this->_fe_u.nbNode ] = this->_wInterp( ig );
            }
        }

        // ALE term: 0.5 div (u^n-w) u v
        mass_divw( - this->density(), _w_loc, _elmatC, this->_fe_u, 0, 0, nc_u );

        // loop on velocity components
        for ( UInt ic = 0;ic < nc_u;ic++ )
        {

             for ( UInt jc = 0;jc < nc_u;jc++ )
             {
                 grad( jc, _elvec, _elmatC, this->_fe_u, this->_fe_u, ic, ic );
                 assemb_mat( _C, _elmatC, this->_fe_u, this->_dof_u, ic, jc );
             }

            // mass
            assemb_mat( _M_u, _elmatM_u, this->_fe_u, this->_dof_u, ic, ic );

            // computing  - (p, \div v) term: the minus sign is in the inner computation
            grad( ic, 1.0, _elmatDtr, this->_fe_u, this->_fe_p, ic, 0 );

            // assembling p div v term and transposed
            assemb_mat_mixed( _trD, _elmatDtr, this->_fe_u, this->_fe_p, this->_dof_u, this->_dof_p, ic, 0 );
            assemb_tr_mat_mixed( 1.0, _D, _elmatDtr, this->_fe_p, this->_fe_u, this->_dof_p, this->_dof_u, 0, ic );
        }
    }

    _trDAux = _trD;
    _CAux = _C;

    // H diag matrix: H= diag( _M_u )/sum( diag( _M_u ) ) where _M_u = mass * rho / dt
    _H = _M_u.giveDiag();
    Real sum = accumulate( _H.begin(), _H.end(), 0.0 );
    // Matrix equilibration
    for ( UInt i = 0; i < _H.size();i++ )
        _H[ i ] = _H[ i ] * dti / sum;

    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;


    // for BC treatment (done at each time-step)
    Real tgv = 1.e02;

    std::cout << "  F-  Applying boundary conditions... ";
    chrono.start();
    _f_u = _f_uWithOutBC;
    this->bcHandler().bdUpdate( this->mesh(), this->_feBd_u, this->_dof_u );
    bcManage( _C, _trD, _f_u, this->mesh(), this->_dof_u, this->bcHandler(), this->_feBd_u, tgv, time );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    //matrices HinvDtr:
    MultInvDiag( _H, _trD, _HinvDtr );
    // ---------------
    // (i) C * V = F_V
    // ---------------
    // AZTEC specifications for each system
    int data_org_i[ AZ_COMM_SIZE ];   // data organisation for C
    int proc_config_i[ AZ_PROC_SIZE ]; // Processor information:
    int options_i[ AZ_OPTIONS_SIZE ]; // Array used to select solver options.
    double params_i[ AZ_PARAMS_SIZE ];   // User selected solver paramters.
    double status_i[ AZ_STATUS_SIZE ];   // Information returned from AZ_solve()
    // indicating success or failure.
    AZ_set_proc_config( proc_config_i, AZ_NOT_MPI );

    //AZTEC matrix and preconditioner
    AZ_MATRIX *C;
    AZ_PRECOND *prec_C;

    int N_eq_i = 3 * this->_dim_u; // number of DOF for each component
    // data_org assigned "by hands" while no parallel computation is performed
    data_org_i[ AZ_N_internal ] = N_eq_i;
    data_org_i[ AZ_N_border ] = 0;
    data_org_i[ AZ_N_external ] = 0;
    data_org_i[ AZ_N_neigh ] = 0;
    data_org_i[ AZ_name ] = DATA_NAME_AZTEC;

    // create matrix and preconditionner
    C = AZ_matrix_create( N_eq_i );
    prec_C = AZ_precond_create( C, AZ_precondition, NULL );

    AZ_set_MSR( C, ( int* ) _pattC.giveRaw_bindx(), ( double* ) _C.giveRaw_value(), data_org_i, 0, NULL, AZ_LOCAL );


    _dataAztec_i.aztecOptionsFromDataFile( options_i, params_i );

    //keep C factorisation and preconditioner reused in my_matvec
    options_i[ AZ_keep_info ] = 1;

    // ---------------
    // (i) C * V = F_V
    // ---------------

    // intermediate velocity computation
    std::cout << "  F-  Solving system (i)... ";
    chrono.start();
    AZ_iterate( this->_u.giveVec(), _f_u.giveVec(), options_i, params_i, status_i,
                proc_config_i, C, prec_C, NULL );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;


    // --------------------------------------------
    // (ii) (D*C^(-1)*trD) * P = D*C^{-1}*F_V = D*V
    // --------------------------------------------
    // AZTEC specifications for the second system
    int proc_config_ii[ AZ_PROC_SIZE ]; // Processor information:
    int options_ii[ AZ_OPTIONS_SIZE ]; // Array used to select solver options.
    double params_ii[ AZ_PARAMS_SIZE ];   // User selected solver paramters.
    double status_ii[ AZ_STATUS_SIZE ];   // Information returned from AZ_solve()
    // indicating success or failure.

    AZ_set_proc_config( proc_config_ii, AZ_NOT_MPI );

    //AZTEC matrix for A_ii=(D*C^{-1}*trD)
    AZ_MATRIX *A_ii;
    AZ_PRECOND *pILU_ii;

    int N_eq_ii = this->_p.size();

    A_ii = AZ_matrix_create( N_eq_ii );
    // data containing the matrices C, D, trD and H as pointers
    // are passed through A_ii and pILU_ii:
    AZ_set_MATFREE( A_ii, &_factor_data,
                    my_matvec < MixedMatr<1, 3, CSRPatt, double>, MixedMatr<3, 1, CSRPatt, double>,
                    std::vector<double>, MSRMatr<double>, Vector > );

    pILU_ii = AZ_precond_create( A_ii, my_precSchur_PC <
                                 MSRMatr<double>,
                                 MixedMatr<1, 3, CSRPatt, double>,
                                 MixedMatr<3, 1, CSRPatt, double>,
                                 std::vector<double>,
                                 MSRMatr<double>,
                                 Vector > , &_factor_data );

    _dataAztec_ii.aztecOptionsFromDataFile( options_ii, params_ii );

    // user preconditioning:
    options_ii[ AZ_precond ] = AZ_user_precond;

    // RHS of the linear system (ii)
    Vector vec_DV( this->_p.size() );


    //matrices HinvC (depends on time):
    MultInvDiag( _H, _C, _HinvC );


    // RHS of the linear system (ii)
    vec_DV = _D * this->_u;

    // case of pure Dirichlet BCs:
    if ( this->bcHandler().hasOnlyEssential())
    {
        vec_DV[ this->_dim_p - 1 ] = 1.0; // correction of the right hand side.
        this->_p[ this->_dim_p - 1 ] = 1.0; // pressure value at the last node.
    }

    std::cout << "  F-  Solving pressure system... ";
    chrono.start();
    AZ_iterate( this->_p.giveVec(), &vec_DV[ 0 ], options_ii, params_ii, status_ii,
                proc_config_ii, A_ii, pILU_ii, NULL );


    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    // ----------------------------
    // (iii) V = V-(C^(-1)*trD) * P
    // ----------------------------

    // everything is done...
    this->_u = this->_u - _invCtrDP;
    std::cout << "  F-  Velocity updated" << std::endl;

    AZ_matrix_destroy( &A_ii );
    AZ_precond_destroy( &pILU_ii );
    AZ_matrix_destroy( &C );
    AZ_precond_destroy( &prec_C );

    _residual_u = _f_uWithOutBC - _CAux * this->_u - _trDAux * this->_p;
}



/*template <typename Mesh>
void NavierStokesAleSolverPC<Mesh>::
iterateTransp( const Real& time )
{

    Chrono chrono;

    _C = _CAux;

    // for BC treatment (done at each time-step)
    Real tgv = 1.e02;

    std::cout << "  F-  Applying boundary conditions... ";
    chrono.start();
    _f_u = _f_uWithOutBC;
    this->BCh_fluid().bdUpdate( this->_mesh, this->_feBd_u, this->_dof_u );
    bcManage( _C, _trD, _f_u, this->_mesh, this->_dof_u, this->BCh_fluid(), this->_feBd_u, tgv, time );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;


    // ---------------
    // (i) C * V = F_V
    // ---------------
    // AZTEC specifications for each system
    int data_org_i[ AZ_COMM_SIZE ];   // data organisation for C
    int proc_config_i[ AZ_PROC_SIZE ]; // Processor information:
    int options_i[ AZ_OPTIONS_SIZE ]; // Array used to select solver options.
    double params_i[ AZ_PARAMS_SIZE ];   // User selected solver paramters.
    double status_i[ AZ_STATUS_SIZE ];   // Information returned from AZ_solve()
    // indicating success or failure.
    AZ_set_proc_config( proc_config_i, AZ_NOT_MPI );

    //AZTEC matrix and preconditioner
    AZ_MATRIX *C;
    AZ_PRECOND *prec_C;

    int N_eq_i = 3 * this->_dim_u; // number of DOF for each component
    // data_org assigned "by hands" while no parallel computation is performed
    data_org_i[ AZ_N_internal ] = N_eq_i;
    data_org_i[ AZ_N_border ] = 0;
    data_org_i[ AZ_N_external ] = 0;
    data_org_i[ AZ_N_neigh ] = 0;
    data_org_i[ AZ_name ] = DATA_NAME_AZTEC;

    // create matrix and preconditionner
    C = AZ_matrix_create( N_eq_i );
    prec_C = AZ_precond_create( C, AZ_precondition, NULL );

    AZ_set_MSR( C, ( int* ) _pattC.giveRaw_bindx(), ( double* ) _C.giveRaw_value(), data_org_i, 0, NULL, AZ_LOCAL );

    _dataAztec_i.aztecOptionsFromDataFile( options_i, params_i );

    //keep C factorisation and preconditioner reused in my_matvec
    options_i[ AZ_keep_info ] = 1;

    // ---------------
    // (i) C * V = F_V
    // ---------------

    // intermediate velocity computation
    std::cout << "  F-  Solving system (i)... ";
    chrono.start();
    AZ_iterate( this->_u.giveVec(), _f_u.giveVec(), options_i, params_i, status_i,
                proc_config_i, C, prec_C, NULL );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;


    // --------------------------------------------
    // (ii) (D*C^(-1)*trD) * P = D*C^{-1}*F_V = D*V
    // --------------------------------------------
    // AZTEC specifications for the second system
    int proc_config_ii[ AZ_PROC_SIZE ]; // Processor information:
    int options_ii[ AZ_OPTIONS_SIZE ]; // Array used to select solver options.
    double params_ii[ AZ_PARAMS_SIZE ];   // User selected solver paramters.
    double status_ii[ AZ_STATUS_SIZE ];   // Information returned from AZ_solve()
    // indicating success or failure.

    AZ_set_proc_config( proc_config_ii, AZ_NOT_MPI );

    //AZTEC matrix for A_ii=(D*C^{-1}*trD)
    AZ_MATRIX *A_ii;
    AZ_PRECOND *pILU_ii;

    int N_eq_ii = _p.size();

    A_ii = AZ_matrix_create( N_eq_ii );
    // data containing the matrices C, D, trD and H as pointers
    // are passed through A_ii and pILU_ii:
    AZ_set_MATFREE( A_ii, &_factor_data,
                    my_matvec < MixedMatr<1, 3, CSRPatt, double>, MixedMatr<3, 1, CSRPatt, double>,
                    std::vector<double>, MSRMatr<double>, Vector > );

    pILU_ii = AZ_precond_create( A_ii, my_precSchur_PC <
                                 MSRMatr<double>,
                                 MixedMatr<1, 3, CSRPatt, double>,
                                 MixedMatr<3, 1, CSRPatt, double>,
                                 std::vector<double>,
                                 MSRMatr<double>,
                                 Vector > , &_factor_data );

    _dataAztec_ii.aztecOptionsFromDataFile( options_ii, params_ii );

    // user preconditioning:
    options_ii[ AZ_precond ] = AZ_user_precond;

    // RHS of the linear system (ii)
    Vector vec_DV( _p.size() );


    // RHS of the linear system (ii)
    vec_DV = _D * _u;

    // case of pure Dirichlet BCs:
    if ( this->BCh_fluid().hasOnlyEssential()
       )
    {
        vec_DV[ this->_dim_p - 1 ] = 1.0; // correction of the right hand side.
        _p[ this->_dim_p - 1 ] = 1.0; // pressure value at the last node.
    }

    std::cout << "  F-  Solving pressure system... ";
    chrono.start();
    AZ_iterate( _p.giveVec(), &vec_DV[ 0 ], options_ii, params_ii, status_ii,
                proc_config_ii, A_ii, pILU_ii, NULL );


    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    // ----------------------------
    // (iii) V = V-(C^(-1)*trD) * P
    // ----------------------------

    // everything is done...
    _u = _u - _invCtrDP;
    std::cout << "  F-  Velocity updated" << std::endl;

    AZ_matrix_destroy( &A_ii );
    AZ_precond_destroy( &pILU_ii );
    AZ_matrix_destroy( &C );
    AZ_precond_destroy( &prec_C );




}*/


//
// Linearized iterate for Newton FSI
//
template <typename Mesh>
void NavierStokesAleSolverPC<Mesh>::
iterateLin( const Real& time, BCHandler& BCh_du )
{

    Chrono chrono;

    // Number of velocity components
    UInt nc_u = this->_u.nbcomp(), iloc, ig;

    std::cout << "  F-  LINEARIZED FLUID SYSTEM\n";

    std::cout << "    F-  Updating right hand side... ";

    //
    // RIGHT HAND SIDE FOR THE LINEARIZED ALE SYSTEM
    //
    chrono.start();

    //initialize right hand side
    _f_duWithOutBC = ZeroVector( _f_duWithOutBC.size() );
    _f_p = ZeroVector( _f_p.size() );

    // Loop on elements
    for ( UInt i = 1; i <= this->mesh().numVolumes(); i++ )
    {

        this->_fe_p.update( this->mesh().volumeList( i ) );
        this->_fe_u.updateFirstDerivQuadPt( this->mesh().volumeList( i ) );

        // initialization of elementary vectors
        _elvec_du.zero();
        _elvec_dp.zero();

        for ( UInt k = 0 ; k < ( UInt ) this->_fe_u.nbNode ; k++ )
        {
            iloc = this->_fe_u.patternFirst( k );
            for ( UInt ic = 0; ic < nc_u; ++ic )
            {
                ig = this->_dof_u.localToGlobal( i, iloc + 1 ) - 1 + ic * this->_dim_u;
                _convect[ iloc + ic * this->_fe_u.nbNode ] = _un( ig ) - this->_wInterp( ig );  // u^n - w^k local
                _w_loc.vec( ) [ iloc + ic * this->_fe_u.nbNode ] = this->_wInterp( ig );                // w^k local
                _uk_loc.vec( ) [ iloc + ic * this->_fe_u.nbNode ] = this->_u( ig );                      // u^k local
                _d_loc.vec( ) [ iloc + ic * this->_fe_u.nbNode ] = this->_dInterp( ig );                // d local
                _dw_loc.vec( ) [ iloc + ic * this->_fe_u.nbNode ] = this->_dwInterp( ig );               // dw local
            }
        }



        for ( UInt k = 0 ; k < ( UInt ) this->_fe_p.nbNode ; k++ )
        {
            iloc = this->_fe_p.patternFirst( k );
            ig = this->_dof_p.localToGlobal( i, iloc + 1 ) - 1;
            _pk_loc[ iloc ] = this->_p( ig );  // p^k local
        }

        //
        // Elementary vectors
        //

        //  - \rho ( -\grad w^k:[I\div d - (\grad d)^T] u^k + ( u^n-w^k )^T[I\div d - (\grad d)^T] (\grad u^k)^T , v  )
        source_mass1( -this->density(), _uk_loc, _w_loc, _convect, _d_loc, _elvec_du, this->_fe_u );

        //  + \rho * ( \grad u^k dw, v  )
        source_mass2( this->density(), _uk_loc, _dw_loc, _elvec_du, this->_fe_u );

        //  - ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  )
        source_stress( -1.0, this->viscosity(), _uk_loc, _pk_loc, _d_loc, _elvec_du, this->_fe_u, this->_fe_p );

        // + \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )
        source_stress2( this->viscosity(), _uk_loc, _d_loc, _elvec_du, this->_fe_u );

        //  + ( (\grad u^k):[I\div d - (\grad d)^T] , q  )
        source_press( 1.0, _uk_loc, _d_loc, _elvec_dp, this->_fe_u, this->_fe_p );

        //
        // Assembling
        //

        // assembling presssure right hand side
        assemb_vec( _f_p, _elvec_dp, this->_fe_p, this->_dof_p, 0 );

        // loop on velocity components
        for ( UInt ic = 0; ic < nc_u; ic++ )
        {
            // assembling velocity right hand side
            assemb_vec( _f_duWithOutBC, _elvec_du, this->_fe_u, this->_dof_u, ic );
        }
    }


    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    // for BC treatment (done at each time-step)
    Real tgv = 1.e02;

    std::cout << "    F-  Applying boundary conditions... ";
    chrono.start();
    _C = _CAux;
    _trD = _trDAux;

    _f_u = _f_duWithOutBC;

    BCh_du.bdUpdate( this->mesh(), this->_feBd_u, this->_dof_u );

    bcManage( _C, _trD, _f_u, this->mesh(), this->_dof_u, BCh_du, this->_feBd_u, tgv, time );

    chrono.stop();
    std::cout << " done in " << chrono.diff() << "s." << std::endl;
    std::cout << "  norm_inf (_f_du) after BC        = " << norm_inf( _f_u ) << std::endl;
    std::cout << "  norm_inf ( difference ) after BC = " << norm_inf( _f_duWithOutBC - _f_u ) << std::endl;

    //matrices HinvDtr:
    MultInvDiag( _H, _trD, _HinvDtr );
    // ---------------
    // (i) C * V = F_V
    // ---------------
    // AZTEC specifications for each system
    int data_org_i[ AZ_COMM_SIZE ];   // data organisation for C
    int proc_config_i[ AZ_PROC_SIZE ]; // Processor information:
    int options_i[ AZ_OPTIONS_SIZE ]; // Array used to select solver options.
    double params_i[ AZ_PARAMS_SIZE ];   // User selected solver paramters.
    double status_i[ AZ_STATUS_SIZE ];   // Information returned from AZ_solve()
    // indicating success or failure.

    AZ_set_proc_config( proc_config_i, AZ_NOT_MPI );

    //AZTEC matrix and preconditioner
    AZ_MATRIX *C;
    AZ_PRECOND *prec_C;

    int N_eq_i = 3 * this->_dim_u; // number of DOF for each component
    // data_org assigned "by hands" while no parallel computation is performed
    data_org_i[ AZ_N_internal ] = N_eq_i;
    data_org_i[ AZ_N_border ] = 0;
    data_org_i[ AZ_N_external ] = 0;
    data_org_i[ AZ_N_neigh ] = 0;
    data_org_i[ AZ_name ] = DATA_NAME_AZTEC;

    // create matrix and preconditionner
    C = AZ_matrix_create( N_eq_i );
    prec_C = AZ_precond_create( C, AZ_precondition, NULL );

    AZ_set_MSR( C, ( int* ) _pattC.giveRaw_bindx(), ( double* ) _C.giveRaw_value(), data_org_i, 0, NULL, AZ_LOCAL );

    _dataAztec_i.aztecOptionsFromDataFile( options_i, params_i );

    //keep C factorisation and preconditioner reused in my_matvec
    options_i[ AZ_keep_info ] = 1;


    // ---------------
    // (i) C * V = F_V
    // ---------------
    options_i[ AZ_recursion_level ] = 1;

    _du = ZeroVector( _du.size() );

    // intermediate velocity computation
    std::cout << "  F-  Solving system (i)... ";
    chrono.start();
    AZ_iterate( _du.giveVec(), _f_u.giveVec(), options_i, params_i, status_i,
                proc_config_i, C, prec_C, NULL );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    //options_i[AZ_recursion_level]=0;

    // --------------------------------------------
    // (ii) (D*C^(-1)*trD) * P = D*C^{-1}*F_V = D*V
    // --------------------------------------------
    // AZTEC specifications for the second system
    int proc_config_ii[ AZ_PROC_SIZE ]; // Processor information:
    int options_ii[ AZ_OPTIONS_SIZE ]; // Array used to select solver options.
    double params_ii[ AZ_PARAMS_SIZE ];   // User selected solver paramters.
    double status_ii[ AZ_STATUS_SIZE ];   // Information returned from AZ_solve()
    // indicating success or failure.

    AZ_set_proc_config( proc_config_ii, AZ_NOT_MPI );

    //AZTEC matrix for A_ii=(D*C^{-1}*trD)
    AZ_MATRIX *A_ii;
    AZ_PRECOND *pILU_ii;

    int N_eq_ii = _dp.size();

    A_ii = AZ_matrix_create( N_eq_ii );
    // data containing the matrices C, D, trD and H as pointers
    // are passed through A_ii and pILU_ii:
    AZ_set_MATFREE( A_ii, &_factor_data_jacobian,
                    my_matvec < MixedMatr<1, 3, CSRPatt, double>, MixedMatr<3, 1, CSRPatt, double>,
                    std::vector<double>, MSRMatr<double>, Vector > );

    pILU_ii = AZ_precond_create( A_ii, my_precSchur_PC <
                                 MSRMatr<double>,
                                 MixedMatr<1, 3, CSRPatt, double>,
                                 MixedMatr<3, 1, CSRPatt, double>,
                                 std::vector<double>,
                                 MSRMatr<double>,
                                 Vector > , &_factor_data_jacobian );

    _dataAztec_ii.aztecOptionsFromDataFile( options_ii, params_ii );

    // user preconditioning:
    options_ii[ AZ_precond ] = AZ_user_precond;

    // RHS of the linear system (ii)
    Vector vec_DV( _dp.size() );


    //matrices HinvC (depends on time):
    MultInvDiag( _H, _C, _HinvC );


    // RHS of the linear system (ii)
    vec_DV = _D * _du - _f_p;

    // case of pure Dirichlet BCs:
    if ( BCh_du.hasOnlyEssential()
       )
    {
        vec_DV[ this->_dim_p - 1 ] = 1.0; // correction of the right hand side.
    }

    _dp = ZeroVector( _dp.size() );


    std::cout << "  F-  Solving pressure system... \n";

//    std::cout << "  norm_inf (vec_DV) = " << norm_inf( vec_DV ) << std::endl;
//    std::cout << "  norm_inf (_f_p) = " << norm_inf( _f_p ) << std::endl;
//    std::cout << "  norm_inf (_D*_du ) = " << norm_inf( _D * _du ) << std::endl;

    chrono.start();
    options_ii[ AZ_recursion_level ] = 1;

    AZ_iterate( _dp.giveVec(), &vec_DV[ 0 ], options_ii, params_ii, status_ii,
                proc_config_ii, A_ii, pILU_ii, NULL );

    chrono.stop();

    std::cout << "done in " << chrono.diff() << " s." << std::endl;
    //options_ii[AZ_recursion_level]=0;

    // ----------------------------
    // (iii) V = V-(C^(-1)*trD) * P
    // ----------------------------
    _du = _du - _invCtrDP;
    std::cout << "  F-  Velocity updated" << std::endl;

    AZ_matrix_destroy( &A_ii );
    AZ_precond_destroy( &pILU_ii );
    AZ_matrix_destroy( &C );
    AZ_precond_destroy( &prec_C );

    _residual_u = _f_duWithOutBC - _CAux * _du - _trDAux * _dp;

    std::cout << "  norm_inf (_residual_du ) = " << norm_inf( _residual_u ) << std::endl;
    std::cout << "  norm_inf (_du )          = " << norm_inf( _du ) << std::endl;
}


template <typename Mesh>
void NavierStokesAleSolverPC<Mesh>::
solveJacobian(  const Real& time, BCHandler& BCh_du )
{

    Chrono chrono;

    // Number of velocity components
    UInt nc_u = this->_u.nbcomp(), iloc, ig;

    std::cout << "  F-  LINEARIZED FLUID SYSTEM\n";

    std::cout << "    F-  Updating right hand side... ";

    //
    // RIGHT HAND SIDE FOR THE LINEARIZED ALE SYSTEM
    //
    chrono.start();

    //initialize right hand side
    _f_duWithOutBC = ZeroVector( _f_duWithOutBC.size() );
    _f_p = ZeroVector( _f_p.size() );

    // Loop on elements
    for ( UInt i = 1; i <= this->mesh().numVolumes(); i++ )
    {

        this->_fe_p.update( this->mesh().volumeList( i ) );
        this->_fe_u.updateFirstDerivQuadPt( this->mesh().volumeList( i ) );

        // initialization of elementary vectors
        _elvec_du.zero();
        _elvec_dp.zero();

        for ( UInt k = 0 ; k < ( UInt ) this->_fe_u.nbNode ; k++ )
        {
            iloc = this->_fe_u.patternFirst( k );
            for ( UInt ic = 0; ic < nc_u; ++ic )
            {
                ig = this->_dof_u.localToGlobal( i, iloc + 1 ) - 1 + ic * this->_dim_u;
                _convect[ iloc + ic * this->_fe_u.nbNode ] = _un( ig ) - this->_wInterp( ig );  // u^n - w^k local
                _w_loc.vec( ) [ iloc + ic * this->_fe_u.nbNode ] = this->_wInterp( ig );                // w^k local
                _uk_loc.vec( ) [ iloc + ic * this->_fe_u.nbNode ] = this->_u( ig );                      // u^k local
                _d_loc.vec( ) [ iloc + ic * this->_fe_u.nbNode ] = this->_dInterp( ig );                // d local
                _dw_loc.vec( ) [ iloc + ic * this->_fe_u.nbNode ] = this->_dwInterp( ig );               // dw local
            }
        }

        for ( UInt k = 0 ; k < ( UInt ) this->_fe_p.nbNode ; k++ )
        {
            iloc = this->_fe_p.patternFirst( k );
            ig = this->_dof_p.localToGlobal( i, iloc + 1 ) - 1;
            _pk_loc[ iloc ] = this->_p( ig );  // p^k local
        }

        //
        // Elementary vectors
        //

        //  - \rho ( \grad( u^n-w^k ):[I\div d - (\grad d)^T] u^k + ( u^n-w^k )^T[I\div d - (\grad d)^T] (\grad u^k)^T , v  )
        source_mass1( -this->density(), _uk_loc, _convect, _d_loc, _elvec_du, this->_fe_u );

        //  + \rho * ( \grad u^k dw, v  )
        source_mass2( this->density(), _uk_loc, _dw_loc, _elvec_du, this->_fe_u );

        //  - ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  )
        source_stress( -1.0, this->viscosity(), _uk_loc, _pk_loc, _d_loc, _elvec_du, this->_fe_u, this->_fe_p );

        // + \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )
        source_stress2( this->viscosity(), _uk_loc, _d_loc, _elvec_du, this->_fe_u );

        //  + ( (\grad u^k):[I\div d - (\grad d)^T] , q  )
        source_press( 1.0, _uk_loc, _d_loc, _elvec_dp, this->_fe_u, this->_fe_p );

        //
        // Assembling
        //

        // assembling presssure right hand side
        assemb_vec( _f_p, _elvec_dp, this->_fe_p, this->_dof_p, 0 );

        // loop on velocity components
        for ( UInt ic = 0; ic < nc_u; ic++ )
        {
            // assembling velocity right hand side
            assemb_vec( _f_duWithOutBC, _elvec_du, this->_fe_u, this->_dof_u, ic );
        }
    }


    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    // for BC treatment (done at each time-step)
    Real tgv = 1.e02;

    std::cout << "    F-  Applying boundary conditions... ";
    chrono.start();
    _C = _CAux;
    _trD = _trDAux;

    _f_u = _f_duWithOutBC;



    BCh_du.bdUpdate( this->mesh(), this->_feBd_u, this->_dof_u );
    bcManage( _C, _trD, _f_u, this->mesh(), this->_dof_u, BCh_du, this->_feBd_u, tgv, time );


    chrono.stop();
    std::cout << " done in " << chrono.diff() << "s." << std::endl;
    std::cout << "  norm_inf (_f_du) after BC= " << norm_inf( _f_u ) << std::endl;
    std::cout << "  norm_inf ( difference ) after BC= " << norm_inf( _f_duWithOutBC - _f_u ) << std::endl;

    //matrices HinvDtr:
    MultInvDiag( _H, _trD, _HinvDtr );
    // ---------------
    // (i) C * V = F_V
    // ---------------
    // AZTEC specifications for each system
    int data_org_i[ AZ_COMM_SIZE ];   // data organisation for C
    int proc_config_i[ AZ_PROC_SIZE ]; // Processor information:
    int options_i[ AZ_OPTIONS_SIZE ]; // Array used to select solver options.
    double params_i[ AZ_PARAMS_SIZE ];   // User selected solver paramters.
    double status_i[ AZ_STATUS_SIZE ];   // Information returned from AZ_solve()
    // indicating success or failure.

    AZ_set_proc_config( proc_config_i, AZ_NOT_MPI );

    //AZTEC matrix and preconditioner
    AZ_MATRIX *C;
    AZ_PRECOND *prec_C;

    int N_eq_i = 3 * this->_dim_u; // number of DOF for each component
    // data_org assigned "by hands" while no parallel computation is performed
    data_org_i[ AZ_N_internal ] = N_eq_i;
    data_org_i[ AZ_N_border ] = 0;
    data_org_i[ AZ_N_external ] = 0;
    data_org_i[ AZ_N_neigh ] = 0;
    data_org_i[ AZ_name ] = DATA_NAME_AZTEC;

    // create matrix and preconditionner
    C = AZ_matrix_create( N_eq_i );
    prec_C = AZ_precond_create( C, AZ_precondition, NULL );

    AZ_set_MSR( C, ( int* ) _pattC.giveRaw_bindx(), ( double* ) _C.giveRaw_value(), data_org_i, 0, NULL, AZ_LOCAL );

    _dataAztec_i.aztecOptionsFromDataFile( options_i, params_i );

    //keep C factorisation and preconditioner reused in my_matvec
    options_i[ AZ_keep_info ] = 1;


    // ---------------
    // (i) C * V = F_V
    // ---------------
    options_i[ AZ_recursion_level ] = 1;

    _du = ZeroVector( _du.size() );

    // intermediate velocity computation
    std::cout << "  F-  Solving system (i)... ";
    chrono.start();
    AZ_iterate( _du.giveVec(), _f_u.giveVec(), options_i, params_i, status_i,
                proc_config_i, C, prec_C, NULL );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    //options_i[AZ_recursion_level]=0;

    // --------------------------------------------
    // (ii) (D*C^(-1)*trD) * P = D*C^{-1}*F_V = D*V
    // --------------------------------------------
    // AZTEC specifications for the second system
    int proc_config_ii[ AZ_PROC_SIZE ]; // Processor information:
    int options_ii[ AZ_OPTIONS_SIZE ]; // Array used to select solver options.
    double params_ii[ AZ_PARAMS_SIZE ];   // User selected solver paramters.
    double status_ii[ AZ_STATUS_SIZE ];   // Information returned from AZ_solve()
    // indicating success or failure.

    AZ_set_proc_config( proc_config_ii, AZ_NOT_MPI );

    //AZTEC matrix for A_ii=(D*C^{-1}*trD)
    AZ_MATRIX *A_ii;
    AZ_PRECOND *pILU_ii;

    int N_eq_ii = _dp.size();

    A_ii = AZ_matrix_create( N_eq_ii );
    // data containing the matrices C, D, trD and H as pointers
    // are passed through A_ii and pILU_ii:
    AZ_set_MATFREE( A_ii, &_factor_data_jacobian,
                    my_matvec < MixedMatr<1, 3, CSRPatt, double>, MixedMatr<3, 1, CSRPatt, double>,
                    std::vector<double>, MSRMatr<double>, Vector > );

    pILU_ii = AZ_precond_create( A_ii, my_precSchur_PC <
                                 MSRMatr<double>,
                                 MixedMatr<1, 3, CSRPatt, double>,
                                 MixedMatr<3, 1, CSRPatt, double>,
                                 std::vector<double>,
                                 MSRMatr<double>,
                                 Vector > , &_factor_data_jacobian );

    _dataAztec_ii.aztecOptionsFromDataFile( options_ii, params_ii );

    // user preconditioning:
    options_ii[ AZ_precond ] = AZ_user_precond;

    // RHS of the linear system (ii)
    Vector vec_DV( _dp.size() );


    //matrices HinvC (depends on time):
    MultInvDiag( _H, _C, _HinvC );


    // RHS of the linear system (ii)
    vec_DV = _D * _du - _f_p;

    // case of pure Dirichlet BCs:
    if ( BCh_du.hasOnlyEssential()
       )
    {
        vec_DV[ this->_dim_p - 1 ] = 1.0; // correction of the right hand side.
    }

    _dp = ZeroVector( _dp.size() );


    std::cout << "  F-  Solving pressure system... \n";
    std::cout << "  norm_inf (vec_DV) = " << norm_inf( vec_DV ) << std::endl;
    std::cout << "  norm_inf (_f_p) = " << norm_inf( _f_p ) << std::endl;
    std::cout << "  norm_inf (_D*_du ) = " << norm_inf( _D * _du ) << std::endl;

    chrono.start();
    options_ii[ AZ_recursion_level ] = 1;

    AZ_iterate( _dp.giveVec(), &vec_DV[ 0 ], options_ii, params_ii, status_ii,
                proc_config_ii, A_ii, pILU_ii, NULL );

    chrono.stop();

    std::cout << "done in " << chrono.diff() << " s." << std::endl;
    //options_ii[AZ_recursion_level]=0;

    // ----------------------------
    // (iii) V = V-(C^(-1)*trD) * P
    // ----------------------------
    _du = _du - _invCtrDP;
    std::cout << "  F-  Velocity updated" << std::endl;

    AZ_matrix_destroy( &A_ii );
    AZ_precond_destroy( &pILU_ii );
    AZ_matrix_destroy( &C );
    AZ_precond_destroy( &prec_C );

    _residual_u = _f_duWithOutBC - _CAux * _du - _trDAux * _dp;

    std::cout << "  norm_inf (_residual_du ) = " << norm_inf( _residual_u ) << std::endl;
    std::cout << "  norm_inf (_du )          = " << norm_inf( _du ) << std::endl;
}




template <typename Mesh>
Vector& NavierStokesAleSolverPC<Mesh>::
residual()
{
    return _residual_u;
}

}
#endif
