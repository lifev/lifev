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
    typedef boost::signal<void ( uint, mesh_type const&, Vector const&, Vector const&, value_type )> iteration_signal_type;

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
        _M_lambda( 0 ),
        _M_source( default_source_term() )
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
        _M_lambda( __nswf._M_lambda ),
        _M_source( __nswf._M_source )
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

    void setSourceTerm( source_type __source )
        {
            _M_source = __source;
        }

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

    void timeAdvance( const Function __source, const Real& __time );

    void iterate( const Real& time );
    //@}



protected:

private:


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

    source_type _M_source;
};

#if 0
Real u0o(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  switch(i) {
  case 1:
  case 2:
    return 0.0;
    break;
  case 3:
     return 0.0;
     break;
  }
  return 0;
}
#endif

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
    _M_solver->timeAdvance( f, time );
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
NavierStokesWithFlux<NSSolver>::timeAdvance( const Function __source, const Real& __time )
{
    // update right hand side for NS solves for new time step
    _M_solver->timeAdvance( __source,__time );

}

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::iterate( const Real& time )
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

    _M_iteration_finish_signal( ( uint )time*100,
                                _M_solver->mesh(),
                                _M_solver->u(),
                                _M_solver->p(),
                                1 );
     //compute the flux of NS: the definitive one
     //
     Qn=_M_solver->flux(_M_fluxes.begin()->first);
     std::cout << "imposed flux" << " " << Q << "\n";
     std::cout << "numerical flux" << " " << Qn << "\n";


}
}
#endif /* __NavierStokesWithFlux_H */
