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
   \author Christian Vergara <>
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

    typedef boost::function<double ( double, double, double )> flux_type;
    typedef std::map<std::string, flux_type> flux_map_type;
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
        _M_ndof( _M_solver->uDof().numTotalDof() )
        {
            Debug( 6010 ) << "constructor from a NavierStokes solver\n";
        }

    NavierStokesWithFlux( NavierStokesWithFlux const & __nswf )
        :
        _M_solver( __nswf._M_solver ),
        _M_ndof( __nswf._M_ndof )
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

    //@}

    /** @name  Mutators
     */
    //@{

    void setFlux( std::string const& __id, flux_type const& __flux )
        {
            Debug( 6010 ) << "imposing flux on boundary" << __id << "\n";
            _M_fluxes[__id]=__flux;
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
    void solve();

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
};

template<typename NSSolver>
void
NavierStokesWithFlux<NSSolver>::solve()
{
    // \warning TO BE REMOVED
    // christian: put your algorithm here the idea is to have one NS
    // solver to avoid duplicating the matrices and you change the BC
    // conditions with modifyBC() which you may have to change if it
    // is not satisfactory, for example:
    // _M_solver->bcHandler().modifyBC( bcname, bcvec);
    // the idea is to pass the id of the boundaries where fluxes will be imposed
    // using the setFlux() memeber function

    // Initialization
    //
    Real dt = _M_solver->timestep();
    Real startT = _M_solver->inittime();
    Real T  = _M_solver->endtime();

    _M_solver->initialize(u0,p0,0.0,dt);

    // simple NS solve that you have to replace with your stuff
    for (Real time=startT+dt ; time <= T; time+=dt)
    {
        _M_solver->timeAdvance(f,time);
        _M_solver->iterate(time);

        _M_iteration_finish_signal( ( uint )time*100,
                                    _M_solver->mesh(),
                                    _M_solver->u(),
                                    _M_solver->p(),
                                    1 );
    }
}
}
#endif /* __NavierStokesWithFlux_H */
