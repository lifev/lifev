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
#include "bcManage.hpp"
#include "bcHandler.hpp"
#include "chrono.hpp"
#include "SolverAztec.hpp"
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
    typedef typename ElasticStructureHandler<Mesh>::source_type source_type;

    //! Constructor
    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature rule
      \param BCh boundary conditions for the displacement
    */
    VenantKirchhofSolver( const GetPot& data_file, const RefFE& refFE, const QuadRule& Qr,
                          const QuadRule& bdQr, BCHandler& BCh );

    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    void timeAdvance( source_type const& source, const Real& time );

    //! Solve the non-linear system
    void iterate();
    void iterate(Vector &_sol);

    //! Output
    void showMe( std::ostream& c = std::cout ) const;

    //! getters

    //! BCHandler getter and setter

    BCHandler const & BC_solid() const {return _BCh;}
    void setBCSolid(const BCHandler & BCd) {_BCh = BCd;}
    //! residual getter
    Vector& residual() {return _residual_d;}

    //! recur setter

    void setRecur(UInt recur) {_recur = recur;}

    void updateJac( Vector& sol, int iter );

    //! solves the tangent problem for newton iterations
    void solveJac( Vector &step, const Vector& res, double& linear_rel_tol);
    //    void solveJac( const Vector& res, double& linear_rel_tol, Vector &step);
    //! solves the tangent problem with custom BC
    void solveJac( Vector &step, const Vector& res, double& linear_rel_tol, BCHandler &BCd );
    void solveLin( Vector &step, const Vector& res, double linear_rel_tol, BCHandler &BCd );
    //! evaluates residual for newton interations
    void evalResidual( Vector &res, const Vector& sol, int iter);

    void setSourceTerm( source_type const& __s ) { _M_source = __s; }
    source_type const& sourceTerm() const { return _M_source; }

    Vector& rhsWithoutBC() { return _rhsWithoutBC; }



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

    //! updates the tangent matrix for newton iterations

    //! files for lists of iterations and residuals per timestep
    std::ofstream _out_iter;
    std::ofstream _out_res;

    //! the present time
    Real _time;

    //! level of recursion for Aztec (has a sens with FSI coupling)
    UInt _recur;

    source_type _M_source;

    //! data for solving tangent problem with aztec
    SolverAztec _linearSolver;
};


//
//                                         IMPLEMENTATION
//
template <typename Mesh>
VenantKirchhofSolver<Mesh>::
VenantKirchhofSolver( const GetPot& data_file, const RefFE& refFE, const QuadRule& Qr,
                      const QuadRule& bdQr, BCHandler& BCh ) :
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
	_out_iter( "out_iter_solid" ),
    _out_res( "out_res_solid" ),
    _time( 0.0 ),
    _recur( 0 )
{

    std::cout << std::endl;
    std::cout << "O-  Displacement unknowns: " << this->_dim << std::endl;
    std::cout << "O-  Computing mass and linear strain matrices... ";

    _linearSolver.setOptionsFromGetPot( data_file, "solid/aztec" );
    _linearSolver.setMatrix( _J );

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
timeAdvance( source_type const& source, const Real& time )
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
    _rhsWithoutBC = ZeroVector( _rhsWithoutBC.size() );

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
    Vector __z = this->_d + this->_dt * _w;
    _rhsWithoutBC += _M * __z;
    _rhsWithoutBC -= _K * this->_d;

    _rhs_w = ( 2.0 / this->_dt ) * this->_d + _w;
    std::cout << std::endl;
    std::cout << "rhsWithoutBC norm = " << norm_2(_rhsWithoutBC) << std::endl;
    std::cout << "_rhs_w norm       = " << norm_2(_rhs_w) << std::endl;
    std::cout << "    _w norm       = " << norm_2(_w) << std::endl;
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

    status = newton( this->_d, *this, norm_inf_adaptor(), _abstol, _reltol, maxiter, _etamax, ( int ) _linesearch, _out_res, _time );

    if ( status == 1 )
    {
        std::ostringstream __ex;
        __ex << "VenantKirchhofSolver::iterate() Inners newton iterations failed to converge\n";
        throw std::logic_error( __ex.str() );
    }
    else
    {
        std::cout << "Number of inner iterations       : " << maxiter << std::endl;
        _out_iter << _time << " " << maxiter << std::endl;
    }

    _w = ( 2.0 / this->_dt ) * this->_d - _rhs_w;

//    evalResidual(_residual_d, _d, 0);
//    _residual_d = -1*(_C*this->_d);
//    std::cout << "rhsWithoutBC norm = " << norm_2(_rhsWithoutBC) << std::endl;
    _residual_d = _C*this->_d - _rhsWithoutBC;
//    _residual_d = -1.*_residual_d;
}

template <typename Mesh>
void VenantKirchhofSolver<Mesh>::
iterate(Vector &_sol)
{

    int status;

    int maxiter = _maxiter;

    status = newton( _sol, *this, norm_inf_adaptor(), _abstol, _reltol, maxiter, _etamax, ( int ) _linesearch, _out_res, _time );

    if ( status == 1 )
    {
        std::ostringstream __ex;
        __ex << "VenantKirchhofSolver::iterate( Vector ) Inners newton iterations failed to converge\n";
        throw std::logic_error( __ex.str() );
    }
    else
    {
        std::cout << "Number of inner iterations       : " << maxiter << std::endl;
        _out_iter << _time << " " << maxiter << std::endl;
    }

    _w = ( 2.0 / this->_dt ) * this->_d - _rhs_w;

    std::cout << "sol norm = " << norm(_sol) << std::endl;

    _residual_d = _C*_sol - _rhsWithoutBC;
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
evalResidual( Vector &res, const Vector& sol, int iter)
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

    bcManageMatrix( _K, _mesh, this->_dof, _BCh, _feBd, 1.0 );

    _rhs = _rhsWithoutBC;

    bcManageVector( _rhs, _mesh, this->_dof, _BCh, _feBd, _time, 1.0 );

    res = _K * sol - _rhs;

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
}



template <typename Mesh>
void VenantKirchhofSolver<Mesh>::
updateJac( Vector& sol, int iter )
{
    std::cout << "    o-  Solid: Updating JACOBIAN in iter " << iter << "  ... ";

    Chrono chrono;
    chrono.start();


    // copy of the linear part
    _J = _C;

    if ( _maxiter > 1 )
    {
        std::cout << " ******** non linear part" << std::endl;

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
//     if (iter == 1)
//     {
//         _J.spy("Jacobian");
//     }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
}




template <typename Mesh>
void VenantKirchhofSolver<Mesh>::
//solveJac( const Vector& res, double& linear_rel_tol, Vector &step)
solveJac( Vector &step, const Vector& res, double& linear_rel_tol)
{
    Chrono chrono;

    _f = res;

    // for BC treatment (done at each time-step)
    Real tgv = 1.0;
    std::cout << "   o-  Applying boundary conditions... ";
    chrono.start();

    // BC manage for the velocity
    if ( !_BCh.bdUpdateDone() )
        _BCh.bdUpdate( _mesh, _feBd, this->_dof );

    bcManageMatrix( _J, _mesh, this->_dof, _BCh, _feBd, tgv );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    _linearSolver.setRecursionLevel( _recur );

    std::cout << "   o-  Solving system... "<< std::flush;
    chrono.start();
    _linearSolver.solve( step , _f);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    //--options[AZ_recursion_level];

//    AZ_matrix_destroy( &J );
//    AZ_precond_destroy( &prec_J );

    _residual_d = _C*step;// - _rhsWithoutBC;
}




template <typename Mesh>
void VenantKirchhofSolver<Mesh>::
solveJac(Vector &step, const Vector& res, double& linear_rel_tol, BCHandler &BCd )
{

    Chrono chrono;

    _f = res;

    // for BC treatment (done at each time-step)
    Real tgv = 1.0;
    std::cout << "  o-  Applying boundary conditions... ";
    chrono.start();

    // BC manage for the velocity
    if ( !BCd.bdUpdateDone() )
        BCd.bdUpdate( _mesh, _feBd, this->_dof );

    bcManageMatrix( _J, _mesh, this->_dof, BCd, _feBd, tgv );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    _linearSolver.setRecursionLevel( _recur );

    std::cout << "  o-  Solving system... "<< std::flush;
    chrono.start();
    _linearSolver.solve( step , _f );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    //--options[AZ_recursion_level];

//    AZ_matrix_destroy( &J );
//    AZ_precond_destroy( &prec_J );

    _residual_d = _C*step;
}


template <typename Mesh>
void VenantKirchhofSolver<Mesh>::
solveLin( Vector &step, const Vector& res, double linear_rel_tol, BCHandler &BCd )
{
    Chrono chrono;

    _f = res;
    _J = _C;

    // for BC treatment (done at each time-step)
    Real tgv = 1.0;
    std::cout << "  S-  Applying boundary conditions... ";
    chrono.start();

    // BC manage for the velocity
    if ( !BCd.bdUpdateDone() )
        BCd.bdUpdate( _mesh, _feBd, this->_dof );

    bcManageMatrix( _J, _mesh, this->_dof, BCd, _feBd, tgv );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    _linearSolver.setRecursionLevel( _recur );

    std::cout << "  o-  Solving system... "<< std::flush;
    chrono.start();
    _linearSolver.solve( step , _f );
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    _w = ( 2.0 / this->_dt ) * step - _rhs_w;

    _residual_d = _C*step;
}



}
#endif
