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

#include "ElasticStructureHandler.hpp"
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "elemOper.hpp"
#include "values.hpp"
#include "pattern.hpp"
#include "assemb.hpp"
#include "bc_manage.hpp"
#include "bcCond.hpp"
#include "chrono.hpp"
#include "dataAztec.hpp"
#include "dataNewton.hpp"
#include "newton.hpp"

namespace LifeV
{

/*!
  \class VenantKirchhofSolver
 
  \brief
  This class solves the linear elastodynamics equations for a (only linear right now)
  St. Venant-Kirchoff material
 
 
*/
template <typename Mesh>
class VenantKirchhofSolver:
            public ElasticStructureHandler<Mesh>, public DataNewton
{

public:

    typedef typename ElasticStructureHandler<Mesh>::Function Function;

    //! Constructor
    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature rule
      \param BCh boundary conditions for the displacement
    */
    VenantKirchhofSolver( const GetPot& data_file, const RefFE& refFE, const QuadRule& Qr,
                          const QuadRule& bdQr, BC_Handler& BCh );

    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    void timeAdvance( const Function source, const Real& time );

    //! Solve the non-linear system
    void iterate();

    //! Output
    void showMe( std::ostream& c = std::cout ) const;

    //! residual getter

    Vector& residual()
    {
        return _residual_d;
    };

    //! friends classes related to the newton solver
    template <class Fct, class Vector, class Real, class Norm>
    friend int newton( Vector& sol, Fct& f, Norm& norm, Real abstol, Real reltol, int& maxit,
                       Real eta_max, int linesearch, std::ofstream& out_res, const Real& time );
    template <class Fct, class Vector, class Real, class Norm>
    friend void lineSearch_cubic( Fct& f, Norm& norm, Vector& residual, Vector& sol, Vector& step,
                                  Real& normRes, Real& lambda, Real slope, int iter );
    template <class Fct, class Vector, class Real, class Norm>
    friend void lineSearch_parab( Fct& f, Norm& norm, Vector& residual, Vector& sol, Vector& step, Real& normRes,
                                  Real& lambda, int iter );

private:

    //! Block pattern of M
    MSRPatt _pattM_block;

    //! Pattern for M
    MixedPattern<3, 3, MSRPatt> _pattM;

    //! Block pattern of K
    MSRPatt _pattK;

    //! Matrix M: mass
    MixedMatr<3, 3, MSRPatt, double> _M;

    MSRMatr<double> _Kl;

    //! Matrix Knl: stiffness non-linear
    MSRMatr<double> _K;

    //! Matrix C: mass + linear stiffness
    MSRMatr<double> _C;

    //! Matrix J: jacobian
    MSRMatr<double> _J;

    //! Elementary matrices and vectors
    ElemMat _elmatK; // stiffnes
    ElemMat _elmatM; // mass
    ElemMat _elmatC; // mass + stiffness
    ElemVec _elvec;  // Elementary right hand side
    ElemVec _dk_loc; // Local displacement

    //! right  hand  side displacement
    PhysVectUnknown<Vector> _rhs;

    //! right  hand  side velocity
    PhysVectUnknown<Vector> _rhs_w;


    //! right  hand  side
    PhysVectUnknown<Vector> _rhsWithoutBC;

    //! right  hand  side
    PhysVectUnknown<Vector> _f;

    //! residual
    PhysVectUnknown<Vector> _residual_d;

    //! data for solving tangent problem with aztec
    DataAztec _dataAztec;

    //! evaluates residual for newton interations
    void evalResidual( Vector&res, const Vector& sol, int iter );

    //! updates the tangent matrix for newton iterations
    void updateJac( Vector& sol, int iter );

    //! solves the tangent problem for newton iterations
    void solveJac( Vector& step, const Vector& res, double& linear_rel_tol );

    //! files for lists of iterations and residuals per timestep
    std::ofstream _out_iter;
    std::ofstream _out_res;

    //! the present time
    Real _time;

    //! level of recursion for Aztec (has a sens with FSI coupling)
    UInt _recur;

    friend class operFS;

};


//
//                                         IMPLEMENTATION
//
template <typename Mesh>
VenantKirchhofSolver<Mesh>::
VenantKirchhofSolver( const GetPot& data_file, const RefFE& refFE, const QuadRule& Qr,
                      const QuadRule& bdQr, BC_Handler& BCh ) :
        ElasticStructureHandler<Mesh>( data_file, refFE, Qr, bdQr, BCh ),
        DataNewton( data_file, "solid/newton" ),
        _pattM_block( this->_dof ),
        _pattM( _pattM_block, "diag" ),
        _pattK( this->_dof, 3 ),
        _M( _pattM ),
        _Kl( _pattK ),
        _K( _pattK ),
        _C( _pattK ),
        _J( _pattK ),
        _elmatK( this->_fe.nbNode, nDimensions, nDimensions ),
        _elmatM( this->_fe.nbNode, nDimensions, nDimensions ),
        _elmatC( this->_fe.nbNode, nDimensions, nDimensions ),
        _elvec( this->_fe.nbNode, nDimensions ),
        _dk_loc( this->_fe.nbNode, nDimensions ),
        _rhs( this->_dim ),
        _rhs_w( this->_dim ),
        _rhsWithoutBC( this->_dim ),
        _f( this->_dim ),
        _residual_d( this->_dim ),
        _dataAztec( data_file, "solid/aztec" ),
        _out_iter( "out_iter_solid" ),
        _out_res( "out_res_solid" ),
        _time( 0.0 ),
        _recur( 0 )
{

    std::cout << std::endl;
    std::cout << "O-  Displacement unknowns: " << this->_dim << std::endl;
    std::cout << "O-  Computing mass and linear strain matrices... ";

    Chrono chrono;
    chrono.start();

    // Matrices initialization
    _M.zeros();
    _Kl.zeros();
    _C.zeros();
    // Number of displacement components
    UInt nc = this->_d.nbcomp();

    //inverse of dt:
    Real dti2 = 2.0 / ( this->_dt * this->_dt );

    // Elementary computation and matrix assembling
    // Loop on elements
    for ( UInt i = 1; i <= _mesh.numVolumes(); i++ )
    {

        this->_fe.updateFirstDerivQuadPt( _mesh.volumeList( i ) );

        _elmatK.zero();
        _elmatM.zero();

        // stiffness
        stiff_strain( _mu, _elmatK, this->_fe );
        stiff_div ( 0.5 * _lambda, _elmatK, this->_fe );

        _elmatC.mat() = _elmatK.mat();

        // mass
        mass( dti2 * _rho, _elmatM, this->_fe, 0, 0, nDimensions );

        _elmatC.mat() += _elmatM.mat();

        // assembling
        for ( UInt ic = 0;ic < nc;ic++ )
        {
            for ( UInt jc = 0;jc < nc;jc++ )
            {
                assemb_mat( _Kl, _elmatK, this->_fe, this->_dof, ic, jc );
                assemb_mat( _C, _elmatC, this->_fe, this->_dof, ic, jc );
            }

            //mass
            assemb_mat( _M, _elmatM, this->_fe, this->_dof, ic, ic );
        }
    }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

}

template <typename Mesh>
void VenantKirchhofSolver<Mesh>::
timeAdvance( const Function source, const Real& time )
{

    UInt ig;

    _time = time;

    std::cout << std::endl;
    std::cout << "O== SOLID: Now we are at time " << _time << " s." << std::endl;

    std::cout << "  o-  Updating mass term on right hand side... ";

    Chrono chrono;
    chrono.start();

    _K = _Kl;

    // Number of displacement components
    UInt nc = this->_d.nbcomp();

    if ( _maxiter > 1 )
    {

        // l`oop on volumes: assembling source term
        for ( UInt i = 1; i <= _mesh.numVolumes(); ++i )
        {

            this->_fe.updateFirstDerivQuadPt( _mesh.volumeList( i ) );

            _elmatK.zero();

            // _dk_loc contains the displacement in the nodes
            for ( UInt j = 0 ; j < ( UInt ) this->_fe.nbNode ; ++j )
            {
                for ( UInt ic = 0; ic < nc; ++ic )
                {
                    ig = this->_dof.localToGlobal( i, j + 1 ) - 1 + ic * this->_dim;
                    _dk_loc[ j + ic * this->_fe.nbNode ] = this->_d( ig );
                }
            }

            // stiffness for non-linear terms
            // 1/2 * \mu * ( [\grad d^k]^T \grad d : \grad v  )
            stiff_dergradbis( _mu * 0.5, _dk_loc, _elmatK, this->_fe );

            // 1/4 * \lambda * ( \tr { [\grad d^k]^T \grad d }, \div v  )
            stiff_derdiv( _lambda * 0.25, _dk_loc, _elmatK, this->_fe );

            for ( UInt ic = 0; ic < nc; ++ic )
            {
                for ( UInt jc = 0;jc < nc;jc++ )
                    assemb_mat( _K, _elmatK, this->_fe, this->_dof, ic, jc );
            }
        }
    }

    // Right hand side for the velocity at time
    _rhsWithoutBC = 0.;

    // loop on volumes: assembling source term
    for ( UInt i = 1; i <= _mesh.numVolumes(); ++i )
    {

        this->_fe.updateFirstDerivQuadPt( _mesh.volumeList( i ) );

        _elvec.zero();

        for ( UInt ic = 0; ic < nc; ++ic )
        {
            compute_vec( source, _elvec, this->_fe, _time, ic ); // compute local vector
            assemb_vec( _rhsWithoutBC, _elvec, this->_fe, this->_dof, ic ); // assemble local vector into global one
        }
    }

    // right hand side without boundary load terms
    _rhsWithoutBC += _M * ( this->_d + this->_dt * _w );
    _rhsWithoutBC -= _K * this->_d;

    _rhs_w = ( 2.0 / this->_dt ) * this->_d + _w;

    //
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

}


template <typename Mesh>
void VenantKirchhofSolver<Mesh>::
iterate()
{

    int status;

    int maxiter = _maxiter;

    status = newton( this->_d, *this, maxnorm, _abstol, _reltol, maxiter, _etamax, ( int ) _linesearch, _out_res, _time );

    if ( status == 1 )
    {
        std::cout << "Inners iterations failed\n";
        exit( 1 );
    }
    else
    {
        std::cout << "Number of inner iterations       : " << maxiter << std::endl;
        _out_iter << _time << " " << maxiter << std::endl;
    }

    _w = ( 2.0 / this->_dt ) * this->_d - _rhs_w;

    _residual_d = _K * this->_d - _rhs;

}



template <typename Mesh>
void VenantKirchhofSolver<Mesh>::
showMe( std::ostream& c ) const
{
    DataElasticStructure<Mesh>::showMe( c );
    c << "\n*** Values for data [solid/newton]\n\n";
    DataNewton::showMe( c );
}

template <typename Mesh>
void VenantKirchhofSolver<Mesh>::
evalResidual( Vector&res, const Vector& sol, int iter )
{


    std::cout << "O-    Computing residual... ";


    Chrono chrono;
    chrono.start();

    // Matrices initialization
    _K = _C;

    if ( _maxiter > 1 )
    {

        UInt ig;

        // Number of displacement components
        UInt nc = this->_d.nbcomp();

        // Elementary computation and matrix assembling
        // Loop on elements
        for ( UInt i = 1; i <= _mesh.numVolumes(); i++ )
        {

            this->_fe.updateFirstDerivQuadPt( _mesh.volumeList( i ) );

            _elmatK.zero();

            // _dk_loc contains the displacement in the nodes
            for ( UInt j = 0 ; j < ( UInt ) this->_fe.nbNode ; ++j )
            {
                for ( UInt ic = 0; ic < nc; ++ic )
                {
                    ig = this->_dof.localToGlobal( i, j + 1 ) - 1 + ic * this->_dim;
                    _dk_loc[ j + ic * this->_fe.nbNode ] = sol( ig );
                }
            }
            // stiffness for non-linear terms

            // 1/2 * \mu * ( [\grad d^k]^T \grad d : \grad v  )
            stiff_dergradbis( _mu * 0.5, _dk_loc, _elmatK, this->_fe );

            // 1/4 * \lambda * ( \tr { [\grad d^k]^T \grad d }, \div v  )
            stiff_derdiv( _lambda * 0.25, _dk_loc , _elmatK, this->_fe );

            // assembling
            for ( UInt ic = 0;ic < nc;ic++ )
                for ( UInt jc = 0;jc < nc;jc++ )
                    assemb_mat( _K, _elmatK, this->_fe, this->_dof, ic, jc );
        }
    }




    if ( !_BCh.bdUpdateDone() )
        _BCh.bdUpdate( _mesh, _feBd, this->_dof );
    bc_manage_matrix( _K, _mesh, this->_dof, _BCh, _feBd, 1.0 );

    _rhs = _rhsWithoutBC;
    bc_manage_vector( _rhs, _mesh, this->_dof, _BCh, _feBd, _time, 1.0 );

    res = _K * sol - _rhs;

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

}



template <typename Mesh>
void VenantKirchhofSolver<Mesh>::
updateJac( Vector& sol, int iter )
{


    std::cout << "    o-  Updating JACOBIAN in iter " << iter << std::endl << "  ... ";

    Chrono chrono;
    chrono.start();


    // copy of the linear part
    _J = _C;

    if ( _maxiter > 1 )
    {

        UInt ig;

        // Number of displacement components
        UInt nc = this->_d.nbcomp();

        // loop on volumes: assembling source term
        for ( UInt i = 1; i <= _mesh.numVolumes(); ++i )
        {

            _fe.updateFirstDerivQuadPt( _mesh.volumeList( i ) );

            _elmatK.zero();

            // _dk_loc contains the displacement in the nodes
            for ( UInt j = 0 ; j < ( UInt ) _fe.nbNode ; ++j )
            {
                for ( UInt ic = 0; ic < nc; ++ic )
                {
                    ig = this->_dof.localToGlobal( i, j + 1 ) - 1 + ic * this->_dim;
                    _dk_loc[ j + ic * _fe.nbNode ] = sol[ ig ];
                }
            }

            // stiffness for non-linear terms
            // 1/2 * \mu * ( [\grad \delta d]^T \grad d^k + [\grad d^k]^T \grad \delta d : \grad v  )
            stiff_dergrad( _mu * 0.5, _dk_loc, _elmatK, _fe );

            // 1/2 * \lambda * ( \tr { [\grad u^k]^T \grad u }, \div v  )
            stiff_derdiv( 0.5 * _lambda, _dk_loc, _elmatK, _fe );

            // assembleing
            for ( UInt ic = 0; ic < nc; ++ic )
                for ( UInt jc = 0; jc < nc; jc++ )
                    assemb_mat( _J, _elmatK, _fe, this->_dof, ic, jc );

        }

    }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

}




template <typename Mesh>
void VenantKirchhofSolver<Mesh>::
solveJac( Vector& step, const Vector& res, double& linear_rel_tol )
{

    Chrono chrono;


    _f = res;

    // for BC treatment (done at each time-step)
    Real tgv = 1.0;
    std::cout << "  o-  Applying boundary conditions... ";
    chrono.start();

    // BC manage for the velocity
    if ( !_BCh.bdUpdateDone() )
        _BCh.bdUpdate( _mesh, _feBd, this->_dof );

    bc_manage_matrix( _J, _mesh, this->_dof, _BCh, _feBd, tgv );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    // AZTEC specifications for the first system
    int data_org[ AZ_COMM_SIZE ];   // data organisation for C
    int proc_config[ AZ_PROC_SIZE ]; // Processor information:
    int options[ AZ_OPTIONS_SIZE ]; // Array used to select solver options.
    double params[ AZ_PARAMS_SIZE ];   // User selected solver paramters.
    double status[ AZ_STATUS_SIZE ];   // Information returned from AZ_solve()
    // indicating success or failure.
    AZ_set_proc_config( proc_config, AZ_NOT_MPI );

    //AZTEC matrix and preconditioner
    AZ_MATRIX *J;
    AZ_PRECOND *prec_J;

    int N_eq = 3 * this->_dim; // number of DOF for each component
    // data_org assigned "by hands" while no parallel computation is performed
    data_org[ AZ_N_internal ] = N_eq;
    data_org[ AZ_N_border ] = 0;
    data_org[ AZ_N_external ] = 0;
    data_org[ AZ_N_neigh ] = 0;
    data_org[ AZ_name ] = 0;

    // create matrix and preconditionner
    J = AZ_matrix_create( N_eq );
    prec_J = AZ_precond_create( J, AZ_precondition, NULL );

    AZ_set_MSR( J, ( int* ) _pattK.giveRaw_bindx(), ( double* ) _J.giveRaw_value(), data_org, 0, NULL, AZ_LOCAL );

    _dataAztec.aztecOptionsFromDataFile( options, params );

    options[ AZ_recursion_level ] = _recur;


    //keep  factorisation and preconditioner reused in my_matvec
    // options_i[AZ_keep_info]= 1;

    //params[AZ_tol]       = linear_rel_tol;

    std::cout << "  o-  Solving system...  ";
    chrono.start();
    AZ_iterate( &step[ 0 ], _f.giveVec(), options, params, status,
                proc_config, J, prec_J, NULL );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    //--options[AZ_recursion_level];

    AZ_matrix_destroy( &J );
    AZ_precond_destroy( &prec_J );

}
}
#endif
