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
   \file NavierStokesWithTwoFluxes.hpp
   \author Christian Vergara <vergara@mate.polimi.it>
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-08
 */

#ifndef __NavierStokesWithTwoFluxes_H
#define __NavierStokesWithTwoFluxes_H 1

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/signal.hpp>

#include <debug.hpp>


namespace LifeV
{
/*!
  \class NavierStokesWithTwoFluxes
  \brief Functor class to impose flux at boundary conditions for Navier-Stokes.

  @author Christian Vergara
  @author Christophe Prud'homme
*/
template<typename NSSolver>
class NavierStokesWithTwoFluxes
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

    NavierStokesWithTwoFluxes( solver_type& __s )
        :
        _M_solver( __s ),
        _M_ndof( _M_solver->uDof().numTotalDof() )
        {
            Debug( 6010 ) << "constructor from a NavierStokes solver\n";
        }

    NavierStokesWithTwoFluxes( NavierStokesWithTwoFluxes const & __nswf )
        :
        _M_solver( __nswf._M_solver ),
        _M_ndof( __nswf._M_ndof )
        {
            Debug( 6010 ) << "copy constructor\n";
        }

    ~NavierStokesWithTwoFluxes()
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
    void solve(int, int);

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
NavierStokesWithTwoFluxes<NSSolver>::solve(int label0, int label1)
{
    // Definition of GMRes alghoritm variables and constants
    //
    Real Q[2],Qn[2];
    Real lambda[2],v1[2],v2[2],r0[2],z[2];
    Real w1[2],w2[2],h[2][2];
    Real absr0;
    PhysVectUnknown<Vector> uzero(_M_solver->uDof().numTotalDof());
    PhysVectUnknown<Vector> uprec(_M_solver->uDof().numTotalDof());
    PhysVectUnknown<Vector> uns1(_M_solver->uDof().numTotalDof());
    PhysVectUnknown<Vector> unso1(_M_solver->uDof().numTotalDof());
    PhysVectUnknown<Vector> unso2(_M_solver->uDof().numTotalDof());
    ScalUnknown<Vector> pns1(_M_solver->pDof().numTotalDof());
    ScalUnknown<Vector> pnso1(_M_solver->pDof().numTotalDof());
    ScalUnknown<Vector> pnso2(_M_solver->pDof().numTotalDof());

    // Initialization
    //
    Real dt = _M_solver->timestep();
    Real startT = _M_solver->inittime();
    Real T  = _M_solver->endtime();

    for (UInt i=0 ; i<_M_solver->uDof().numTotalDof(); ++i) {
        uzero[i]=0.;
        uprec[i]=0.;
    }
    // Change the BC for NS1
    //
    UInt dim_lambda = _M_solver->uDof().numTotalDof();
    Vector vec_lambda0(dim_lambda);
    Vector vec_lambda1(dim_lambda);
    BCVector bcvec0(vec_lambda0,dim_lambda,1);
    BCVector bcvec1(vec_lambda1,dim_lambda,1);
    _M_solver->bcHandler().modifyBC(label0,bcvec0);
    _M_solver->bcHandler().modifyBC(label1,bcvec1);

    // Initialize
    //
    _M_solver->initialize(u0,p0,0.0,dt);

    // Temporal loop
    //
    for (Real time=startT+dt ; time <= T; time+=dt)
    {
        Q[0]=(_M_fluxes[label0])(time);
        Q[1]=(_M_fluxes[label1])(time);
        lambda[0]=-Q[0];
        lambda[1]=-Q[1];
        std::cout << "imposed flux 0" << " " << Q[0] << "\n";
        std::cout << "imposed flux 1" << " " << Q[1] << "\n";

        // NS1
        //
        // Costruction of vectors vec_lambda
        //
        for (UInt i=0 ; i < dim_lambda; ++i) {
            vec_lambda0[i]=-lambda[0];
            vec_lambda1[i]=lambda[1];
        }

        //std::cout << "\n";
        std::cout << "start NS1 time" << " " << time << "\n";
        _M_solver->timeAdvance(f,time);
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
        std::cout << "flux 0 before the update" << " " << Qn[0] << "\n";
        std::cout << "flux 1 before the update" << " " << Qn[1] << "\n";

        // GMRes algorithm
        //
        r0[0]=Qn[0]-Q[0];
        r0[1]=Qn[1]-Q[1];
        absr0=std::sqrt(r0[0]*r0[0]+r0[1]*r0[1]);
        v1[0]=r0[0]/absr0;
        v1[1]=r0[1]/absr0;

        // NSo1
        //
        // Changing BC
        //
        for (UInt i=0 ; i < dim_lambda; ++i) {
            vec_lambda0[i]=v1[0];
            vec_lambda1[i]=v1[1];
        }

        std::cout << "start NSo1 time" << " " << time << "\n";

        _M_solver->u()=uzero;
        _M_solver->timeAdvance(f,time);
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

        std::cout << "coeff" << " " << w1[0] << "\n";
        std::cout << "coeff" << " " << w1[1] << "\n";

        h[0][0]=w1[0]*v1[0]+w1[1]*v1[1];
        w1[0]=w1[0]-h[0][0]*v1[0];
        w1[1]=w1[1]-h[0][0]*v1[1];
        h[1][0]=std::sqrt(w1[0]*w1[0]+w1[1]*w1[1]);
        //h[1,0]=w1[0]*w1[0]+w1[1]*w1[1];
        v2[0]=w1[0]/h[1][0];
        v2[1]=w1[1]/h[1][0];

        std::cout << "coeff" << " " << h[0][0] << "\n";
        std::cout << "coeff" << " " << h[1][0] << "\n";
        std::cout << "coeff" << " " << w1[0] << "\n";
        std::cout << "coeff" << " " << w1[1] << "\n";
        std::cout << "coeff" << " " << v1[0] << "\n";
        std::cout << "coeff" << " " << v1[1] << "\n";
        std::cout << "coeff" << " " << v2[0] << "\n";
        std::cout << "coeff" << " " << v2[1] << "\n";


        // NSo2
        //
        // Changing BC
        //
        for (UInt i=0 ; i < dim_lambda; ++i) {
            vec_lambda0[i]=v2[0];
            vec_lambda1[i]=v2[1];
        }

        std::cout << "start NSo2 time" << " " << time << "\n";

        _M_solver->u()=uzero;
        _M_solver->timeAdvance(f,time);
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

        h[0][1]=w2[0]*v1[0]+w2[1]*v1[1];
        w2[0]=w2[0]-h[0][1]*v1[0];
        w2[1]=w2[1]-h[0][1]*v1[1];
        h[1][1]=w2[0]*v2[0]+w2[1]*v2[1];
        w2[0]=w2[0]-h[1][1]*v2[0];
        w2[1]=w2[1]-h[1][1]*v2[1];

        // check
        //
        std::cout << "w2" << " " << w2[0] << "\n";
        std::cout << "w2" << " " << w2[1] << "\n";
        std::cout << "coeff" << " " << h[0][0] << "\n";
        std::cout << "coeff" << " " << h[1][0] << "\n";
        std::cout << "coeff" << " " << h[0][1] << "\n";
        std::cout << "coeff" << " " << h[1][1] << "\n";

        // Update lambda
        //
        z[0]=-h[1][1]*absr0/(h[1][0]*h[0][1]-h[0][0]*h[1][1]);
        z[1]=-h[1][0]*absr0/(h[1][0]*h[0][1]-h[0][0]*h[1][1]);
        lambda[0]=lambda[0]+v1[0]*z[0]+v2[0]*z[1];
        lambda[1]=lambda[1]+v1[1]*z[0]+v2[1]*z[1];

        std::cout << "coeff" << " " << z[0] << "\n";
        std::cout << "coeff" << " " << z[1] << "\n";

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
    }
}
}
#endif /* __NavierStokesWithTwoFluxes_H */
