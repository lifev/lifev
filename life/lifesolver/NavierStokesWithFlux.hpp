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
NavierStokesWithFlux<NSSolver>::solve(int label)
{
    // Definition of GMRes alghoritm variables and constants
    //
    Real lambda0,Q,Qno,Qn;
    Real lambda;
    Real r0,y,z;
    Real absr0;
    int v;
   
    // Initialization
    //
    Real dt = _M_solver->timestep();
    Real startT = _M_solver->inittime();
    Real T  = _M_solver->endtime();

    // Stationary Navier-Stokes (NSo)
    //
    _M_solver->initialize(u0o,p0,0.0,dt);

    for (Real time=startT+dt ; time <= startT+dt; time+=dt) {
      std::cout << "\n";
      std::cout << "start NSo\n";
      _M_solver->timeAdvance(f,time);
      _M_solver->iterate(time);
    }

    // Store the solutions of NSo
    //
    PhysVectUnknown<Vector> u_nso(_M_solver->uDof().numTotalDof());
    PhysVectUnknown<Vector> p_nso(_M_solver->pDof().numTotalDof());
    u_nso=_M_solver->u();
    for (UInt jj=0 ; jj<_M_solver->pDof().numTotalDof(); jj+=1) {
      p_nso[jj]=_M_solver->p()[jj];
    }

    // compute the flux of NSo  
    //  
    Qno=_M_solver->flux(label); 
    std::cout << "\n";
    std::cout << "end NSo\n";
    
    // Change the BC for the non stationary NS
    //
    UInt dim_lambda = _M_solver->uDof().numTotalDof();
    Vector vec_lambda(dim_lambda);
    BCVector bcvec(vec_lambda,dim_lambda,1);
    _M_solver->bcHandler().modifyBC(label,bcvec);

    // Navier Stokes in temporal loop
    //
    _M_solver->initialize(u0,p0,0.0,dt);
    
    // Temporal loop
    //
    for (Real time=startT+dt ; time <= T; time+=dt)
    {
      Q=(_M_fluxes[label])(time);
      lambda0=-Q;

      // Costruction of vector vec_lambda
      //
      for (UInt i=0 ; i < dim_lambda; ++i) {
        vec_lambda[i]=-lambda0;
      }
  
      std::cout << "\n";
      std::cout << "start NS time" << " " << time << "\n";
      _M_solver->timeAdvance(f,time);
      _M_solver->iterate(time);

      _M_iteration_finish_signal( ( uint )time*100,
                                    _M_solver->mesh(),
                                    _M_solver->u(),
                                    _M_solver->p(),
                               1 );

      //compute the flux of NS
      //
      Qn=_M_solver->flux(label); 
      std::cout << "flux before update" << " " << Qn << "\n"; 

      // compute the variables to update lambda
      r0=Qn-Q;
      if(r0>0){
        absr0=r0;
        v=1;
      }
      else{
        absr0=-r0;
        v=-1;
      }
      y=absr0/Qno;
      z=v*y;
      lambda=lambda0+z;

      // update the velocity and the pressure
      _M_solver->u()=_M_solver->u()-z*u_nso;
      for (UInt jj=0 ; jj<_M_solver->pDof().numTotalDof(); jj+=1) {
        _M_solver->p()[jj]=_M_solver->p()[jj]-z*p_nso[jj];
      }

      // Save the final solutions
      // 
      _M_solver->postProcess();   

      //compute the flux of NS: the definitive one
      //
      Qn=_M_solver->flux(label);
      std::cout << "imposed flux" << " " << Q << "\n"; 
      std::cout << "numerical flux" << " " << Qn << "\n";
    }
}
}
#endif /* __NavierStokesWithFlux_H */
