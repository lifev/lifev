/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-12

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
   \file main_ns_flux.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-12
 */
#include <life/lifecore/life.hpp>
#include <life/lifesolver/NavierStokesSolverPC.hpp>
#include <life/lifesolver/NavierStokesWithFlux.hpp>
#include <life/lifecore/chrono.hpp>
#include <ud_functions.hpp>
#include <life/lifecore/GetPot.hpp>

namespace LifeV
{
struct toEnsight
{
    void
    operator()( Real __time,
                RegionMesh3D<LinearTetra> const& __mesh,
                Vector const& __u,
                Vector const& __p,
                Real __flux ) const
        {
            std::cout << "============================================================\n"
                      << "Call ensight writer here to save the data \n"
                      << "============================================================\n";

        }
};
}

int
main(int argc, char** argv)
{
    using namespace LifeV;
    using namespace std;

    // Reading from data file
    //
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);

    // BC common to all the NS
    //
    BCFunctionBase u_wall(u1);
    BCFunctionBase out_flow(u1);
    //BCFunctionBase in_flow(u1); // needs for two fluxes imposed at outlets
    BCFunctionBase in_flow(uo);
    BCHandler BCh_u(5);
    BCh_u.addBC("Wall",   2, Essential, Full, u_wall,  3);
    BCh_u.addBC("InFlow", 1, Natural,   Full, in_flow, 3);
    BCh_u.addBC("OutFlow", 3, Natural,   Full, out_flow, 3);
    BCh_u.addBC("InFlowWall", 4, Essential,   Full, out_flow, 3);
    BCh_u.addBC("OutFlowWall", 5, Essential,   Full, out_flow, 3);

    // Navier-Stokes Solver
    //
    typedef NavierStokesSolverPC< RegionMesh3D<LinearTetra> > ns_type;
    boost::shared_ptr<ns_type> __ns ( new ns_type(data_file,
                                                  feTetraP1bubble, feTetraP1,
                                                  quadRuleTetra15pt,quadRuleTria3pt,
                                                  quadRuleTetra5pt, quadRuleTria3pt,
                                                  BCh_u) );
    __ns->showMe();
    __ns->setSourceTerm( f );
    __ns->yesFlux(1);

    NavierStokesWithFlux<ns_type> __ns_with_flux( __ns );

    // Impose the fluxes for initialize
    //
    __ns_with_flux.setFlux(1, my_flux_cost);
    //__ns_with_flux.setFlux(3, my_flux_cos2);

    //Set the strategy: 0 for the exact version; 1 for the inexact one. If flux imposed is one, the two versions are the same
    //
    __ns_with_flux.setStrategy(0);

    toEnsight EnsightFilter;
    __ns_with_flux.doOnIterationFinish( EnsightFilter  );

    Real dt = __ns_with_flux.timestep();
    Real startT = __ns_with_flux.inittime();
    Real T  = __ns_with_flux.endtime();

    __ns_with_flux.initialize(u0,p0,0.0,dt);

    //ofstream outfile("flusso_inesatto.txt");

    for (Real time=startT+dt ; time <= T; time+=dt)
    {

       // Impose the fluxes
       //
       __ns_with_flux.setFlux(1, my_flux_cost); //costant
       //__ns_with_flux.setFlux(3, my_flux_cos2); //costant
       //__ns_with_flux.setFlux(1, my_flux_cos); //cosinusoidal
       //__ns_with_flux.setFlux(1, my_flux_physio); // physiological

       __ns_with_flux.iterate( time );
       }

    return EXIT_SUCCESS;
}

