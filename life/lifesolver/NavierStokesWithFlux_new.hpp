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

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/signal.hpp>

#include <debug.hpp>


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
    //typedef boost::function< Real ( Real, Real, Real Real, ID ) > source_type;
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
        _M_u_nso( _M_solver->uDof().numTotalDof() ),
        _M_p_nso( _M_solver->pDof().numTotalDof() ),
        _M_Qno( 0 ),
        _M_vec_lambda( _M_solver->uDof().numTotalDof() ),
        _M_lambda( 0 )
        {
            Debug( 6010 ) << "constructor from a NavierStokes solver\n";
        }

    NavierStokesWithFlux( NavierStokesWithFlux const & __nswf )
        :
        _M_solver( __nswf._M_solver ),
        _M_ndof( __nswf._M_ndof ),
        _M_u_nso( __nswf._M_u_nso ),
        _M_p_nso( __nswf._M_p_nso ),
        _M_Qno( __nswf._M_Qno ),
        _M_vec_lambda( __nswf._M_vec_lambda ),
        _M_lambda( __nswf._M_lambda )
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

    flux_type flux( std::string const& __id ) const
        {
            flux_map_const_iterator __it = _M_fluxes.find( __id );
            if (  __it == _M_fluxes.end() )
            {
                std::ostringstream __ex;
                __ex << "Invalid bounday identifier for flux imposition: "
                     << __id;
                throw std::invalid_argument( __ex.str() );
            }
            return *__it;
        }

    Real timestep() const { return _M_solver->timestep();}
    Real inittime() const { return _M_solver->inittime();}
    Real endtime() const { return _M_solver->endtime();}

    Real pressure() const
        {
            return _M_lambda;
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
    void solve(int);

    void initialize( const Function& u0, const Function& p0, Real t0, Real dt );

    void timeAdvance( source_type const& __source, const Real& __time );

    void iterate( const Real& time );
    //@}



protected:

private:

    //! one flux version of the algorithm
    void iterate_one_flux( Real const& );

    //! two fluxes version of the algorithm
    void iterate_two_fluxes( Real const& );

private:

    //! Navier-Stokes solver
    solver_type _M_solver;

    //! number of degree of freedom for Lagrange multipliers
    UInt _M_ndof;

    //! flux map
    flux_map_type _M_fluxes;

    //! signal at each iteration
    iteration_signal_type _M_iteration_finish_signal;

    //! Velocity correction for Flux imposition
    PhysVectUnknown<Vector> _M_u_nso;

    //! Pressure correction for Flux imposition
    ScalUnknown<Vector> _M_p_nso;

    //! Flux from NSo
    Real _M_Qno;

    //! lagrange multiplier to impose flux
    Vector _M_vec_lambda;

    Real _M_lambda;

};

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::initialize( const Function& u0, const Function& p0, Real t0, Real dt )
{
    Debug( 6010 ) << "start NSo\n";

    // if one flux do the NSo solves
    // Stationary Navier-Stokes (NSo)
    //
    _M_solver->initialize(u0o,p0,0.0,dt);

    Real startT = _M_solver->inittime();
    Real time=startT+dt;
    _M_solver->timeAdvance( _M_solver->sourceTerm(), time );
    _M_solver->iterate( time );

    // Store the solutions of NSo
    //
    _M_u_nso=_M_solver->u();
    _M_p_nso=_M_solver->p();

    // compute the flux of NSo
    //
    _M_Qno=_M_solver->flux(_M_fluxes.begin()->first);
    Debug( 6010 ) << "end NSo\n";

    //
    // Change the BC for the non stationary NS
    //

    BCVector bcvec( _M_vec_lambda, _M_vec_lambda.size(), 1 );
    _M_solver->bcHandler().modifyBC( _M_fluxes.begin()->first, bcvec );

    // Navier Stokes in temporal loop
    //
    _M_solver->initialize(u0,p0,0.0,dt);

}
template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::timeAdvance( source_type const& __source, Real const& __time )
{
    // update right hand side for NS solves for new time step
    _M_solver->timeAdvance( __source,__time );

}

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::iterate( const Real& time )
{
    switch ( _M_fluxes.size() == 1 )
    {
        case 1:
            timeAdvance( _M_solver->sourceTerm(), time );
            iterate_one_flux( time );
        break;
        case 2:
            iterate_two_fluxes( time );
            break;
    }
}
template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::iterate_one_flux( const Real& time )
{
    // update the left hand side and solve the system at time \c time
    Real Q=_M_fluxes.begin()->second( time );
    Real lambda0 = -Q;

    // update vec_lambda
    _M_vec_lambda = ScalarVector( _M_vec_lambda.size(), -lambda0 );

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
    _M_solver->u()-=z * _M_u_nso;
    _M_solver->p()-=z * _M_p_nso;

    _M_iteration_finish_signal( time,
                                _M_solver->mesh(),
                                _M_solver->u(),
                                _M_solver->p(),
                                1 );
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
    Debug( 6010 ) << "starting two fluxes version at time " << time << "\n";
#if 0
    int label0 = _M_fluxes.begin().first;
    int label1 = boost::next( _M_fluxes.begin() ).first;
    Q[0] = _M_fluxes[label0](time);
    Q[1] = _M_fluxes[label1](time);
    lambda[0]=-Q[0];
    lambda[1]=-Q[1];
    std::cout << "imposed flux 0" << " " << Q[0] << "\n";
    std::cout << "imposed flux 1" << " " << Q[1] << "\n";

    // NS1
    //
    // Costruction of vectors vec_lambda
    //
    vec_lambda0 = ScalarVector( dim_lambda, -lambda[0] );
    vec_lambda1 = ScalarVector( dim_lambda,  lambda[1] );

    //std::cout << "\n";
    Debug( 6010 ) << "start NS1 time" << " " << time << "\n";
    _M_solver->timeAdvance(_M_solver->sourceTerm(),time);
    _M_solver->iterate(time);

    //_M_iteration_finish_signal( ( uint )time*100,
    //                             _M_solver->mesh(),
    //                             _M_solver->u(),
    //                             _M_solver->p(),
    //                        1 );



    // Store the solution
    //
    uns1=_M_solver->u();
    pns1=_M_solver->p();



    // Compute the fluxes of NS1
    //
    Qn[0]=_M_solver->flux(label0);
    Qn[1]=_M_solver->flux(label1);
    debug( 6010 ) << "flux 0 before the update" << " " << Qn[0] << "\n";
    Debug( 6010 ) << "flux 1 before the update" << " " << Qn[1] << "\n";

    // GMRes algorithm
    //
    r0 = Qn - Q;
    absr0 = norm_2( r0 );
    v1 = r0/absr0;

    // NSo1
    //
    // Changing BC
    //
    vec_lambda0 = ScalarVector( dim_lambda, v1[0] );
    vec_lambda1 = ScalarVector( dim_lambda, v1[1] );

    Debug( 6010 ) << "start NSo1 time" << " " << time << "\n";

    _M_solver->u()=uzero;
    _M_solver->timeAdvance(_M_solver->sourceTerm(),time);
    //_M_solver->f()=uzero;
    _M_solver->u()=uprec;
    _M_solver->iterate(time);


    // _M_iteration_finish_signal( ( uint )time*100,
    //                             _M_solver->mesh(),
    //                            _M_solver->u(),
    //                            _M_solver->p(),
    //                       1 );

    // Store the solution
    //
    unso1=_M_solver->u();
    pnso1=_M_solver->p();

    // Compute the fluxes of NSo1
    //
    w1[0]=_M_solver->flux(label0);
    w1[1]=_M_solver->flux(label1);

    std::cout << "w1" << " " << w1 << "\n";

    h( 0, 0 ) = sum( w1*v1 );
    w1 -= h( 0, 0 )*v1;
    h( 1, 0 ) = norm_2( w1 );
    v2 = w1/h1( 1, 0 );

    std::cout << "h  = " << h << "\n";
    std::cout << "v1 = " << v1 << "\n";
    std::cout << "v2 = " << v2 << "\n";
    std::cout << "w1 = " << w1 << "\n";


    // NSo2
    //
    // Changing BC
    //
    vec_lambda0 = ScalarVector( dim_lambda, v2[0] );
    vec_lambda1 = ScalarVector( dim_lambda, v2[1] );

    Debug( 6010 ) << "start NSo2 time" << " " << time << "\n";

    _M_solver->u()=uzero;
    _M_solver->timeAdvance(_M_solver->sourceTerm(),time);
    //_M_solver->f()=uzero;
    _M_solver->u()=uprec;
    _M_solver->iterate(time);

    //      _M_iteration_finish_signal( ( uint )time*100,
    //                            _M_solver->mesh(),
    //                            _M_solver->u(),
    //                            _M_solver->p(),
    //                       1 );

    // Store the solution
    //
    unso2=_M_solver->u();
    pnso2=_M_solver->p();

    // Compute the fluxes of NSo2
    //
    w2[0]=_M_solver->flux(label0);
    w2[1]=_M_solver->flux(label1);

    h( 0, 1 ) = sum( w2*v1 );
    w2 -= h( 0, 1 )*v1;
    h( 1, 1 ) = sum( w2*v2 );
    w2 -= h( 1, 1 )*v2;

    //
    // check
    //
    std::cout << "w2 = " << w2 << "\n"
              << "h  = " << h << "\n";

    //
    // Update lambda
    //
    Real __v = absr0/( h( 1, 0 )*h( 0, 1 ) - h( 0, 0 )*h( 1, 1 ) );
    z[0] = -h( 1, 1 )*__v;
    z[1] = -h( 1, 0 )*__v;
    lambda += v1*z[0] + v2*z[1];
    std::cout << "z      = " << z << "\n"
              << "lambda = " << lambda << "\n";

    //
    // Update the velocity and the pressure
    //
    _M_solver->u()=uns1-z[0]*unso1-z[1]*unso2;
    _M_solver->p()=pns1-z[0]*pnso1-z[1]*pnso2;

    // Store the velocity
    //
    uprec=_M_solver->u();

    // Save the final solutions
    //
    _M_solver->postProcess();

    // Compute the fluxes of NS: the definitive one
    //
    Qn[0]=_M_solver->flux(label0);
    Qn[1]=_M_solver->flux(label1);
    std::cout << "imposed flux 0" << " " << Q[0] << "\n";
    std::cout << "numerical flux 0" << " " << Qn[0] << "\n";
    std::cout << "imposed flux 1" << " " << Q[1] << "\n";
    std::cout << "numerical flux 1" << " " << Qn[1] << "\n";
#endif

}
}
#endif /* __NavierStokesWithFlux_H */
