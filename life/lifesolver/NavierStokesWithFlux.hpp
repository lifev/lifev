/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): christian Vergara <>
             Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-08

  Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano

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
   \file NavierStokesWithFlux.hpp
   \author Christian Vergara <vergara@mate.polimi.it>
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-08
 */

#ifndef __NavierStokesWithFlux_H
#define __NavierStokesWithFlux_H 1

#include <stdexcept>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/signal.hpp>
//#include <ensight7Writer.hpp>

#include <life/lifecore/debug.hpp>


namespace LifeV
{
namespace
{
struct default_source_term
{
    Real operator()( Real __t, Real __x,  Real __y,  Real __z, ID __id )
        {
            return 0;
        }
};
struct default_velocity
{
    Real operator()( Real __t, Real __x,  Real __y,  Real __z, ID __id )
        {
            return 0;
        }
};
enum FluxStrategies
  {
    STRATEGY_FLUX_GMRES_EXACT,        /**< exact flux computation albeit require several NS solves */
    STRATEGY_FLUX_GMRES_INEXACT,      /**< inexact flux computation but cheaper */
  };

}
/*!
  \class NavierStokesWithFlux
  \brief Functor class to impose flux at boundary conditions for Navier-Stokes.

  @author Christian Vergara
  @author Christophe Prud'homme
*/
template<typename NSSolver>
class NavierStokesWithFlux
{
public:


    /** @name Typedefs
     */
    //@{

    typedef Real value_type;

    typedef boost::shared_ptr<NSSolver> solver_type;

    typedef boost::function< Real ( Real ) > flux_type;
    typedef boost::function< Real ( Real, Real, Real, Real, ID ) > source_type;
    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );

    typedef std::map<int, flux_type> flux_map_type;
    typedef typename flux_map_type::iterator flux_map_iterator;
    typedef typename flux_map_type::const_iterator flux_map_const_iterator;

    typedef typename NSSolver::mesh_type mesh_type;
    typedef boost::signal<void ( Real,
                                 mesh_type const&,
                                 PhysVectUnknown<Vector> const&, ScalUnknown<Vector> const&,
                                 value_type )> iteration_signal_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    NavierStokesWithFlux( solver_type& __s )
        :

        _M_solver( __s ),
        _M_ndof( _M_solver->uDof().numTotalDof() ),
        _M_unso( _M_solver->uDof().numTotalDof() ),
        _M_uzero(_M_solver->uDof().numTotalDof()),
        _M_uprec(_M_solver->uDof().numTotalDof()),
        _M_uns1(_M_solver->uDof().numTotalDof()),
        _M_unso1(_M_solver->uDof().numTotalDof()),
        _M_unso2(_M_solver->uDof().numTotalDof()),
        _M_unso1_staz(_M_solver->uDof().numTotalDof()),
        _M_unso2_staz(_M_solver->uDof().numTotalDof()),
        _M_pns1(_M_solver->pDof().numTotalDof()),
        _M_pnso1(_M_solver->pDof().numTotalDof()),
        _M_pnso2(_M_solver->pDof().numTotalDof()),
        _M_pnso1_staz(_M_solver->pDof().numTotalDof()),
        _M_pnso2_staz(_M_solver->pDof().numTotalDof()),
        _M_pnso( _M_solver->pDof().numTotalDof() ),
        _M_Qno( 0 ),
        _M_vec_lambda( _M_solver->uDof().numTotalDof() ),
        _M_lambda( 0 ),
        _M_lambda0( 0 ),
        _M_lambda1( 0 )
        {
            Debug( 6010 ) << "constructor from a NavierStokes solver\n";
        }

    NavierStokesWithFlux( NavierStokesWithFlux const & __nswf )
        :
        _M_solver( __nswf._M_solver ),
        _M_ndof( __nswf._M_ndof ),
        _M_unso( __nswf._M_unso ),
        _M_uzero(__nswf._M_uzero),
        _M_uprec(__nswf._M_uprec),
        _M_uns1(__nswf._M_uns1),
        _M_unso1(__nswf._M_unso1),
        _M_unso2(__nswf._M_unso2),
        _M_unso1_staz(__nswf._M_unso1_staz),
        _M_unso2_staz(__nswf._M_unso2_staz),
        _M_pns1(__nswf._M_pns1),
        _M_pnso1(__nswf._M_pnso1),
        _M_pnso2(__nswf._M_pnso2),
        _M_pnso1_staz(__nswf._M_pnso1_staz),
        _M_pnso2_staz(__nswf._M_pnso2_staz),
        _M_pnso( __nswf._M_pnso ),
        _M_Qno( __nswf._M_Qno ),
        _M_vec_lambda( __nswf._M_vec_lambda ),
        _M_lambda( __nswf._M_lambda ),
        _M_lambda0( __nswf._M_lambda0 ),
        _M_lambda1( __nswf._M_lambda1 )
        {
            Debug( 6010 ) << "copy constructor\n";
        }

    ~NavierStokesWithFlux()
        {
            Debug( 6010 ) << "destructor\n";
        }

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    //flux_type flux( std::string const& __id ) const
    flux_type flux( int lab ) const
        {
            flux_map_const_iterator __it = _M_fluxes.find( lab );
            if (  __it == _M_fluxes.end() )
            {
                std::ostringstream __ex;
                __ex << "Invalid bounday identifier for flux imposition: "
                     << lab;
                throw std::invalid_argument( __ex.str() );
            }
            return *__it;
        }

    Real timestep() const { return _M_solver->timestep();}
    Real inittime() const { return _M_solver->inittime();}
    Real endtime() const { return _M_solver->endtime();}

    /** @ mean normal stress on section gamma where is imposed
	the flux Q (iterate_one_flux)
     */
    Real pressure() const
        {
            return _M_lambda;
        }
    //@}

    /** @ mean normal stress on section gamma0 where is imposed
	the flux Q(0) (iterate_two_fluxes)
     */
    Real pressure0() const
        {
            return _M_lambda0;
        }
    //@}

    /** @ mean normal stress on section gamma1 where is imposed
	the flux Q(1) (iterate_two_fluxes)
     */
    Real pressure1() const
        {
            return _M_lambda1;
        }
    //@}

    /** @name  Mutators
     */
    //@{

    void setFlux( int lab, flux_type const& __flux )
        {
            Debug( 6010 ) << "imposing flux on boundary" << lab << "\n";
            _M_fluxes[lab]=__flux;
        }

    //@}

  void setStrategy( FluxStrategies strategy ){
      thestrategy = strategy;
  }

  FluxStrategies strategy() { return thestrategy; }



    /** @name  Methods
     */
    //@{

    template<typename Observer>
    void doOnIterationFinish( Observer& __observer )
        {
            Debug( 6010 ) << "adding observer to iteration finish signal\n";
            _M_iteration_finish_signal.connect( __observer );
        }

    /**
       Use Lagrange multipliers to impose the flux
    */

    void initialize( const Function& u0, const Function& p0, Real t0, Real dt );

    void iterate( const Real& time );

    void iterate( const Real& time, const int& jj );


    //@}



protected:

private:

    //! one flux version of the algorithm
    void iterate_one_flux( Real const& );
    void iterate_one_flux( Real const&, int const& );

    //! two fluxes version of the algorithm
    void iterate_two_fluxes( Real const& );

    //! two fluxes inexact version of the algorithm
    void iterate_two_fluxes_inexact( Real const& );

    //! one flux version of the algorithm
    void initialize_one_flux( const Function& , const Function& , Real , Real  );

    //! two fluxes version of the algorithm
    void initialize_two_fluxes( const Function& , const Function& , Real , Real  );

    //! two fluxes inexact version of the algorithm
    void initialize_two_fluxes_inexact( const Function& , const Function& , Real , Real  );


private:

    //! Navier-Stokes solver
    solver_type _M_solver;

    //! number of degree of freedom for Lagrange multipliers
    UInt _M_ndof;

    //! flux map
    flux_map_type _M_fluxes;

    //! signal at each iteration
    iteration_signal_type _M_iteration_finish_signal;

    //! Variables ror the algorithm
    PhysVectUnknown<Vector> _M_unso;
    PhysVectUnknown<Vector> _M_uzero;
    PhysVectUnknown<Vector> _M_uprec;
    PhysVectUnknown<Vector> _M_uns1;
    PhysVectUnknown<Vector> _M_unso1;
    PhysVectUnknown<Vector> _M_unso2;
    PhysVectUnknown<Vector> _M_unso1_staz;
    PhysVectUnknown<Vector> _M_unso2_staz;
    ScalUnknown<Vector> _M_pns1;
    ScalUnknown<Vector> _M_pnso1;
    ScalUnknown<Vector> _M_pnso2;
    ScalUnknown<Vector> _M_pnso1_staz;
    ScalUnknown<Vector> _M_pnso2_staz;
    ScalUnknown<Vector> _M_pnso;

    //! Flux from NSo
    Real _M_Qno;

    //! lagrange multiplier to impose flux
    Vector _M_vec_lambda;

    Real _M_lambda,_M_lambda0,_M_lambda1;

    FluxStrategies thestrategy;
  //std::ofstream outfile;
};

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::initialize( const Function& u0, const Function& p0, Real t0, Real dt )
{

    switch ( _M_fluxes.size() )
    {
        case 1:
            initialize_one_flux(u0,p0,t0,dt);
 	    break;
        case 2:
      if (thestrategy==STRATEGY_FLUX_GMRES_EXACT){
            initialize_two_fluxes(u0,p0,t0,dt);
      }
          if (thestrategy==STRATEGY_FLUX_GMRES_INEXACT){
            initialize_two_fluxes_inexact(u0,p0,t0,dt);
      }
            break;
        default:
            std::ostringstream __ex;
            __ex << "The number of flux is invalid it is : " << _M_fluxes.size() << "\n"
                 << "you have to specify either one flux or two fluxes for this algorithm to work\n"
                 << "using the setFlux( label, flux ) member function\n";
            throw std::logic_error( __ex.str() );
            break;
    }
}

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::initialize_one_flux( const Function& u0, const Function& p0, Real /*t0*/, Real dt )
{
    Debug( 6020 ) << "start NSo\n";

    // if one flux do the NSo solves
    // Stationary Navier-Stokes (NSo)
    //

    BCVector bcvec( _M_vec_lambda, _M_vec_lambda.size(), 1 );
    _M_solver->bcHandler().modifyBC( _M_fluxes.begin()->first, bcvec );
    _M_vec_lambda = ScalarVector( _M_vec_lambda.size(), 1 );

    _M_solver->initialize(u0o,p0,0.0,dt);

    Real startT = _M_solver->inittime();
    Real time=startT+dt;
    _M_solver->timeAdvance( _M_solver->sourceTerm(), time );
    _M_solver->iterate( time );

    // Store the solutions of NSo
    //
    _M_unso=_M_solver->u();
    _M_pnso=_M_solver->p();

    // compute the flux of NSo
    //
    _M_Qno=_M_solver->flux(_M_fluxes.begin()->first);
    Debug( 6020 ) << "end NSo\n";

    //
    // Change the BC for the non stationary NS
    //

    _M_solver->bcHandler().modifyBC( _M_fluxes.begin()->first, bcvec );

    // Navier Stokes in temporal loop
    //
    _M_solver->initialize(u0,p0,0.0,dt);

}

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::initialize_two_fluxes( const Function& u0, const Function& p0, Real /*t0*/, Real dt )
{
    _M_solver->initialize(u0,p0,0.0,dt);
}

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::initialize_two_fluxes_inexact( const Function& u0, const Function& p0, Real /*t0*/, Real dt )
{
    Debug( 6020 ) << "start NSo1\n";

    int label0 = _M_fluxes.begin()->first;
    int label1 = boost::next( _M_fluxes.begin() )->first;
    UInt dim_lambda=_M_solver->uDof().numTotalDof();
    Vector vec_lambda0(dim_lambda);
    Vector vec_lambda1(dim_lambda);
    BCVector bcvec0(vec_lambda0,dim_lambda,1);
    BCVector bcvec1(vec_lambda1,dim_lambda,1);
    vec_lambda0 = ScalarVector( vec_lambda0.size(),1 );
    vec_lambda1 = ScalarVector( vec_lambda1.size(),0 );

    // Change the BC for the stationary NSo1
    //
    _M_solver->bcHandler().modifyBC(label0,bcvec0 );
    _M_solver->bcHandler().modifyBC(label1,bcvec1 );

    // Stationary Navier-Stokes (NSo1)
    //
    _M_solver->initialize(u0o,p0,0.0,dt);

    Real startT = _M_solver->inittime();
    Real time=startT+dt;
    _M_solver->timeAdvance( _M_solver->sourceTerm(), time );

    _M_solver->iterate( time );


    // Store the solution of NSo1
    //
    _M_unso1_staz=_M_solver->u();
    _M_pnso1_staz=_M_solver->p();

    Debug( 6020 ) << "end NSo1\n";

    // Change the BC for the stationary NSo2
    //
    _M_solver->bcHandler().modifyBC(label0,bcvec1 );
    _M_solver->bcHandler().modifyBC(label1,bcvec0 );

    Debug( 6020 ) << "start NSo2\n";

    // Stationary Navier-Stokes (NSo2)
    //
    _M_solver->initialize(u0o,p0,0.0,dt);

    _M_solver->timeAdvance( _M_solver->sourceTerm(), time );
    _M_solver->iterate( time );

    // Store the solutions of NSo2
    //
    _M_unso2_staz=_M_solver->u();

    Debug( 6020 ) << "end NSo2\n";

    // Initialize the non stationary NS
    //
    _M_solver->initialize(u0,p0,0.0,dt);
}

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::iterate( const Real& time )
{
    switch ( _M_fluxes.size() )
    {
        case 1:
            iterate_one_flux( time );
        break;
        case 2:
          if (thestrategy==STRATEGY_FLUX_GMRES_EXACT){
            iterate_two_fluxes( time );
      }
          if (thestrategy==STRATEGY_FLUX_GMRES_INEXACT){
            iterate_two_fluxes_inexact( time );
      }
            break;
        default:
            std::ostringstream __ex;
            __ex << "The number of flux is invalid it is : " << _M_fluxes.size() << "\n"
                 << "you have to specify either one flux or two fluxes for this algorithm to work\n"
                 << "using the setFlux( label, flux ) member function\n";
            throw std::logic_error( __ex.str() );
            break;
    }
}

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::iterate( const Real& time, const int& jj )
{
    switch ( _M_fluxes.size() )
    {
        case 1:
            iterate_one_flux( time, jj );
        break;
        case 2:
          if (thestrategy==STRATEGY_FLUX_GMRES_EXACT){
            iterate_two_fluxes( time );
      }
          if (thestrategy==STRATEGY_FLUX_GMRES_INEXACT){
            iterate_two_fluxes_inexact( time );
      }
            break;
        default:
            std::ostringstream __ex;
            __ex << "The number of flux is invalid it is : " << _M_fluxes.size() << "\n"
                 << "you have to specify either one flux or two fluxes for this algorithm to work\n"
                 << "using the setFlux( label, flux ) member function\n";
            throw std::logic_error( __ex.str() );
            break;
    }
}

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::iterate_one_flux( const Real& time )
{
    std::cout << "start NS time" << " " << time << "\n";

    // update the left hand side and solve the system at time \c time
    Real Q=_M_fluxes.begin()->second( time );
    Real lambda0 = -Q;

    Debug( 6020 ) << "imposed flux" << " " << Q << "\n";

    // update vec_lambda
    _M_vec_lambda = ScalarVector( _M_vec_lambda.size(), -lambda0 );

    _M_solver->timeAdvance( _M_solver->sourceTerm(), time );
    _M_solver->iterate(time);

    //compute the flux of NS
    //
    Real Qn=_M_solver->flux( _M_fluxes.begin()->first );
    Debug( 6020 ) << "flux before update" << " " << Qn << "\n";

    // compute the variables to update lambda
    Real r0 = Qn-Q;
    Real v;
    Real absr0;
    if(r0>0)
    {
        absr0=r0;
        v=1;
    }
    else
    {
        absr0=-r0;
        v=-1;
    }
    Real y = absr0/_M_Qno;
    Real z = v*y;

    _M_lambda=lambda0+z;

    // update the velocity and the pressure
    _M_solver->u()-=z * _M_unso;
    _M_solver->p()-=z * _M_pnso;

    _M_iteration_finish_signal( time,
                                _M_solver->mesh(),
                                _M_solver->u(),
                                _M_solver->p(),
                                1 );
     // Save the final solutions
     //
    //     _M_solver->postProcess();
    //    outensight7Mesh3D( _M_solver->mesh(), _M_solver->u(), _M_solver->p(), time);

     //compute the flux of NS: the definitive one
     //
     Qn=_M_solver->flux(_M_fluxes.begin()->first);
     Debug( 6020 ) << "imposed flux" << " " << Q << "\n";
     Debug( 6020 ) << "numerical flux" << " " << Qn << "\n";

}

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::iterate_one_flux( const Real& time, const int& jj )
{
    Debug( 6010 ) << "start NS time" << " " << time << "\n";

    // update the left hand side and solve the system at time \c time
    Real Q=_M_fluxes.begin()->second( jj );
    Real lambda0 = -Q;

    Debug( 6010 ) << "imposed flux" << " " << Q << "\n";

    // update vec_lambda
    _M_vec_lambda = ScalarVector( _M_vec_lambda.size(), -lambda0 );

    _M_solver->timeAdvance( _M_solver->sourceTerm(), time );
    _M_solver->iterate(time);

    //compute the flux of NS
    //
    Real Qn=_M_solver->flux( _M_fluxes.begin()->first );
    Debug( 6010 ) << "flux before update" << " " << Qn << "\n";

    // compute the variables to update lambda
    Real r0 = Qn-Q;
    Real v;
    Real absr0;
    if(r0>0)
    {
        absr0=r0;
        v=1;
    }
    else
    {
        absr0=-r0;
        v=-1;
    }
    Real y = absr0/_M_Qno;
    Real z = v*y;

    _M_lambda=lambda0+z;

    // update the velocity and the pressure
    _M_solver->u()-=z * _M_unso;
    _M_solver->p()-=z * _M_pnso;

    _M_iteration_finish_signal( time,
                                _M_solver->mesh(),
                                _M_solver->u(),
                                _M_solver->p(),
                                1 );
     // Save the final solutions
     //
    //     _M_solver->postProcess();
    //    outensight7Mesh3D( _M_solver->mesh(), _M_solver->u(), _M_solver->p(), time);

     //compute the flux of NS: the definitive one
     //
     Qn=_M_solver->flux(_M_fluxes.begin()->first);
     Debug( 6010 ) << "imposed flux" << " " << Q << "\n";
     Debug( 6010 ) << "numerical flux" << " " << Qn << "\n";

}

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::iterate_two_fluxes( const Real& time )
{
    Debug( 6020 ) << "starting two fluxes version at time " << time << "\n";

    int label0 = _M_fluxes.begin()->first;
    int label1 = boost::next( _M_fluxes.begin() )->first;
    Debug( 6020 ) << "imposing fluxes on BC : " << label0 << " and BC : "<< label1 << "\n";

    Vector Q( 2 ),Qn( 2 );
    Vector lambda( 2 ),v1( 2 ),v2( 2 ),r0( 2 ),z( 2 );
    Vector w1( 2 ),w2( 2 );
    Matrix h( 2, 2 );
    Real absr0;


    UInt dim_lambda=_M_solver->uDof().numTotalDof();
    Vector vec_lambda0(dim_lambda);
    Vector vec_lambda1(dim_lambda);
    BCVector bcvec0(vec_lambda0,dim_lambda,1);
    BCVector bcvec1(vec_lambda1,dim_lambda,1);
    _M_solver->bcHandler().modifyBC(label0,bcvec0);
    _M_solver->bcHandler().modifyBC(label1,bcvec1);

    _M_uzero = ZeroVector( _M_uzero.size() );

    Q[0] = _M_fluxes[label0](time);
    Q[1] = _M_fluxes[label1](time);
    lambda = Q;

    _M_lambda0 = lambda[0];
    _M_lambda1 = lambda[1];

    std::cout << "imposed flux 0" << " " << Q[0] << "\n";
    std::cout << "imposed flux 1" << " " << Q[1] << "\n";

    // NS1
    //
    // Construction of vectors vec_lambda
    //
    vec_lambda0 = ScalarVector( vec_lambda0.size(), lambda[0] );
    vec_lambda1 = ScalarVector( vec_lambda1.size(), lambda[1] );

    std::cout << "start NS1 time" << " " << time << "\n";
    _M_solver->timeAdvance(_M_solver->sourceTerm(),time);
    _M_solver->iterate(time);

    _M_iteration_finish_signal( ( uint )time*100,
                                 _M_solver->mesh(),
                                 _M_solver->u(),
                                 _M_solver->p(),
                            1 );

    // Store the solution
    //
    _M_uns1=_M_solver->u();
    _M_pns1=_M_solver->p();

    // Compute the fluxes of NS1
    //
    Qn[0]=_M_solver->flux(label0);
    Qn[1]=_M_solver->flux(label1);
    Debug( 6020 ) << "flux 0 NS1" << " " << Qn[0] << "\n";
    Debug( 6020 ) << "flux 1 NS1" << " " << Qn[1] << "\n";

    // GMRes algorithm
    //
    r0 = Qn - Q;
    absr0 = norm_2( r0 );
    v1 = r0/absr0;

    // NSo1
    //
    // Changing BC
    //
    vec_lambda0 = ScalarVector( vec_lambda0.size(), -v1[0] );
    vec_lambda1 = ScalarVector( vec_lambda1.size(), -v1[1] );

    std::cout << "start NSo1 time" << " " << time << "\n";

    _M_solver->u()=_M_uzero;
    _M_solver->timeAdvance(_M_solver->sourceTerm(),time);
    _M_solver->u()=_M_uprec;
    _M_solver->iterate(time);


    _M_iteration_finish_signal( ( uint )time*100,
                                 _M_solver->mesh(),
                                _M_solver->u(),
                                _M_solver->p(),
                           1 );

    // Store the solution
    //
    _M_unso1=_M_solver->u();
    _M_pnso1=_M_solver->p();

    // Compute the fluxes of NSo1
    //
    w1[0]=_M_solver->flux(label0);
    w1[1]=_M_solver->flux(label1);

    h( 0, 0 ) = inner_prod( w1, v1 );
    w1 -= h( 0, 0 )*v1;
    h( 1, 0 ) = norm_2( w1 );
    v2 = w1/h( 1, 0 );

    // NSo2
    //
    // Changing BC
    //
    vec_lambda0 = ScalarVector( vec_lambda0.size(), -v2[0] );
    vec_lambda1 = ScalarVector( vec_lambda1.size(), -v2[1] );

    std::cout << "start NSo2 time" << " " << time << "\n";

    _M_solver->u()=_M_uzero;
    _M_solver->timeAdvance(_M_solver->sourceTerm(),time);
    _M_solver->u()=_M_uprec;
    _M_solver->iterate(time);

    _M_iteration_finish_signal( ( uint )time*100,
                               _M_solver->mesh(),
                               _M_solver->u(),
                                _M_solver->p(),
                           1 );

    // Store the solution
    //
    _M_unso2=_M_solver->u();
    _M_pnso2=_M_solver->p();

    // Compute the fluxes of NSo2
    //
    w2[0]=_M_solver->flux(label0);
    w2[1]=_M_solver->flux(label1);

    h( 0, 1 ) = inner_prod( w2, v1 );
    w2 -= h( 0, 1 )*v1;
    h( 1, 1 ) = inner_prod( w2, v2 );
    w2 -= h( 1, 1 )*v2;

    // check
    //
    std::cout << "w2 = " << w2 << "\n";

    //
    // Update lambda
    //
    Real tt = absr0/( h( 1, 0 )*h( 0, 1 )-h( 0, 0 )*h( 1, 1 ) );
    z[0] = -h( 1, 1 )*tt;
    z[1] =  h( 1, 0 )*tt;
    lambda += v1*z[0]+v2*z[1];
    std::cout << "lambda = " << lambda << "\n";

    _M_lambda0 = lambda[0];
    _M_lambda1 = lambda[1];


    //
    // Update the velocity and the pressure
    //
    std::cout << "update the velocity and the pressure\n";
    _M_solver->u()=_M_uns1-z[0]*_M_unso1-z[1]*_M_unso2;
    _M_solver->p()=_M_pns1-z[0]*_M_pnso1-z[1]*_M_pnso2;

    // Store the velocity
    //
    _M_uprec=_M_solver->u();

    //    outensight7Mesh3D( _M_solver->mesh(), _M_solver->u(), _M_solver->p(), time);

    // Save the final solutions
    //
    //    _M_solver->postProcess();

    // Compute the fluxes of NS: the definitive one
    //
    Qn[0]=_M_solver->flux(label0);
    Qn[1]=_M_solver->flux(label1);
    std::cout << "imposed flux 0" << " " << Q[0] << "\n";
    std::cout << "numerical flux 0" << " " << Qn[0] << "\n";
    std::cout << "imposed flux 1" << " " << Q[1] << "\n";
    std::cout << "numerical flux 1" << " " << Qn[1] << "\n";
}

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::iterate_two_fluxes_inexact( const Real& time )
{
    Debug( 6020 ) << "starting inexact two fluxes version at time " << time << "\n";

    int label0 = _M_fluxes.begin()->first;
    int label1 = boost::next( _M_fluxes.begin() )->first;
    Debug( 6020 ) << "imposing fluxes on BC : " << label0 << " and BC : "<< label1 << "\n";

    Vector Q( 2 ),Qn( 2 );
    Vector lambda( 2 ),v1( 2 ),v2( 2 ),r0( 2 ),z( 2 );
    Vector w1( 2 ),w2( 2 );
    Matrix h( 2, 2 );
    Real absr0;

    UInt dim_lambda=_M_solver->uDof().numTotalDof();
    Vector vec_lambda0(dim_lambda);
    Vector vec_lambda1(dim_lambda);
    BCVector bcvec0(vec_lambda0,dim_lambda,1);
    BCVector bcvec1(vec_lambda1,dim_lambda,1);
    _M_solver->bcHandler().modifyBC(label0,bcvec0);
    _M_solver->bcHandler().modifyBC(label1,bcvec1);

    _M_uzero = ZeroVector( _M_uzero.size() );

    Q[0] = _M_fluxes[label0](time);
    Q[1] = _M_fluxes[label1](time);
    lambda = Q;

    _M_lambda0 = lambda[0];
    _M_lambda1 = lambda[1];

    std::cout << "imposed flux 0" << " " << Q[0] << "\n";
    std::cout << "imposed flux 1" << " " << Q[1] << "\n";

    // NS
    //
    // Construction of vectors vec_lambda
    //
    vec_lambda0 = ScalarVector( vec_lambda0.size(), lambda[0] );
    vec_lambda1 = ScalarVector( vec_lambda1.size(), lambda[1] );

    std::cout << "start NS time" << " " << time << "\n";
    _M_solver->timeAdvance(_M_solver->sourceTerm(),time);
    _M_solver->iterate(time);

    _M_iteration_finish_signal( ( uint )time*100,
                                 _M_solver->mesh(),
                                 _M_solver->u(),
                                 _M_solver->p(),
                            1 );

    // Store the solution
    //
    _M_uns1=_M_solver->u();
    _M_pns1=_M_solver->p();

    // Compute the fluxes of NS
    //
    Qn[0]=_M_solver->flux(label0);
    Qn[1]=_M_solver->flux(label1);
    Debug( 6020 ) << "flux 0 before update" << " " << Qn[0] << "\n";
    Debug( 6020 ) << "flux 1 before update" << " " << Qn[1] << "\n";

    // GMRes algorithm
    //
    r0 = Qn - Q;
    absr0 = norm_2( r0 );
    v1 = r0/absr0;

    // Compute real (inexact) solution of NSo1
    //
    _M_unso1=v1[0]*_M_unso1_staz+v1[1]*_M_unso2_staz;
    _M_pnso1=v1[0]*_M_pnso1_staz+v1[1]*_M_pnso2_staz;

    // Compute the fluxes of NSo1
    //
    _M_solver->u()=_M_unso1;
    w1[0]=_M_solver->flux(label0);
    w1[1]=_M_solver->flux(label1);

    h( 0, 0 ) = inner_prod( w1, v1 );
    w1 -= h( 0, 0 )*v1;
    h( 1, 0 ) = norm_2( w1 );
    v2 = w1/h( 1, 0 );

    // Compute real (inexact) solution of NSo2
    //
    _M_unso2=v2[0]*_M_unso1_staz+v2[1]*_M_unso2_staz;
    _M_pnso2=v2[0]*_M_pnso1_staz+v2[1]*_M_pnso2_staz;

    // Compute the fluxes of NSo2
    //
    _M_solver->u()=_M_unso2;
    w2[0]=_M_solver->flux(label0);
    w2[1]=_M_solver->flux(label1);

    h( 0, 1 ) = inner_prod( w2, v1 );
    w2 -= h( 0, 1 )*v1;
    h( 1, 1 ) = inner_prod( w2, v2 );
    w2 -= h( 1, 1 )*v2;

    // check
    //
    std::cout << "w2 = " << w2 << "\n";

    // Update lambda
    //
    Real tt = absr0/( h( 1, 0 )*h( 0, 1 )-h( 0, 0 )*h( 1, 1 ) );
    z[0] = -h( 1, 1 )*tt;
    z[1] =  h( 1, 0 )*tt;
    lambda += v1*z[0]+v2*z[1];
    std::cout << "lambda = " << lambda << "\n";

    _M_lambda0 = lambda[0];
    _M_lambda1 = lambda[1];


    // Update the velocity and the pressure
    //
    _M_solver->u()=_M_uns1-z[0]*_M_unso1-z[1]*_M_unso2;
    _M_solver->p()=_M_pns1-z[0]*_M_pnso1-z[1]*_M_pnso2;

    // Save the final solutions
    //
    //    _M_solver->postProcess();

    // Compute the fluxes of NS: the definitive one
    //
    Qn[0]=_M_solver->flux(label0);
    Qn[1]=_M_solver->flux(label1);
    std::cout << "imposed flux 0" << " " << Q[0] << "\n";
    std::cout << "numerical flux 0" << " " << Qn[0] << "\n";
    std::cout << "imposed flux 1" << " " << Q[1] << "\n";
    std::cout << "numerical flux 1" << " " << Qn[1] << "\n";
    //outfile.open("flusso_esatto.txt",std::ios::app);
    //outfile << Q[0] << "\n";
    //outfile << Qn[0] << "\n";
    //outfile << Q[1] << "\n";
    //outfile << Qn[1] << "\n";
    //outfile.close();
}
}
#endif /* __NavierStokesWithFlux_H */
