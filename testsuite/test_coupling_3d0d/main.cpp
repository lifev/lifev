/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include "lifeV.hpp"
#include "NavierStokesSolverPC.hpp"
#include "chrono.hpp"
#include "ud_functions.hpp"
#include "GetPot.hpp"
#include "zeroDModelSolver.hpp"
#include "NavierStokesWithFlux_new.hpp"

#include <ensight7Writer.hpp>
#include <medit_wrtrs.hpp>

namespace LifeV
{
struct toEnsight
{
    void
    operator()( Real __time,
                RegionMesh3D<LinearTetra> const& __mesh,
                PhysVectUnknown<Vector> const& __u,
                ScalUnknown<Vector> const& __p,
                Real __flux ) const
        {
            std::cout << "============================================================\n"
                      << "Call ensight writer here to save the data \n"
                      << "============================================================\n";
            outensight7Mesh3D( __mesh, __u, __p, __time );
        }
};
struct toMedit
{
    toMedit( GetPot __data )
        :
        _M_data( __data )
        {
        }
    void
    operator()( Real __time,
                RegionMesh3D<LinearTetra> const& __mesh,
                PhysVectUnknown<Vector> const& __u,
                ScalUnknown<Vector> const& __p,
                Real __flux ) const
        {
            std::cout << "============================================================\n"
                      << "Call medit writer here to save the data at time: " << __time << "\n"
                      << "============================================================\n";
            std::ostringstream name;
            name <<  __time*1000;
            int __dim_u = __u.size()/__u.nbcomp();
            // postprocess data file for medit
            wr_medit_ascii_scalar( "press." + name.str() + ".bb", __p.giveVec(), __p.size() );
            wr_medit_ascii_scalar( "vel_x." + name.str() + ".bb", __u.giveVec(), __mesh.numVertices() );
            wr_medit_ascii_scalar( "vel_y." + name.str() + ".bb", __u.giveVec() + __dim_u, __mesh.numVertices() );
            wr_medit_ascii_scalar( "vel_z." + name.str() + ".bb", __u.giveVec() + 2 * __dim_u, __mesh.numVertices() );

            std::string __dir = _M_data( "fluid/discretization/mesh_dir", "." );
            std::string __file = _M_data( "fluid/discretization/mesh_file", "mesh.mesh" );

            system( ( "ln -s " + __dir + __file + " press." + name.str() + ".mesh" ).data() );
            system( ( "ln -s " + __dir + __file + " vel_x." + name.str() + ".mesh" ).data() );
            system( ( "ln -s " + __dir + __file + " vel_y." + name.str() + ".mesh" ).data() );
            system( ( "ln -s " + __dir + __file + " vel_z." + name.str() + ".mesh" ).data() );
        }
private:
    GetPot _M_data;
};

class flux_adaptor
{
public:

    flux_adaptor( double __Q )
        :
        _M_Q( __Q )
        {}

    Real operator()( Real /*time*/ )
        {
            return _M_Q;
        }
private:
    Real _M_Q;
};
}
/**
    This test couples the Navier-Stokes equations (3D Model) with a lumped parameter model
    (0D Model - electric network).
    Two couple strategies ares used:

        Mean pressure : the 3D Model gives the flux information to the 0D Model wereas
                        the 0D Model passes the pressure to the 3D.

	Flow rate     : the 3D Model gives the pressure information to the 0D Model wereas
                        the 0D Model passes the flux to the 3D.

*/

int main(int argc, char** argv)
{
    using namespace LifeV;

    // Reading from data file
    //
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);

    int strategy = data_file("problem/strategy", 1 );

    BCFunctionBase u_wall(u1);
    BCFunctionBase out_flow(u1);
    BCHandler BCh_u(5);


    // Navier-Stokes Solver
    //
    typedef NavierStokesSolverPC< RegionMesh3D<LinearTetra> > NS;
    boost::shared_ptr<NS> ns( new NS (data_file, feTetraP1bubble, feTetraP1,
                                      quadRuleTetra15pt,quadRuleTria3pt,
                                      quadRuleTetra5pt, quadRuleTria3pt, BCh_u) );
    ns->showMe();
    ns->setSourceTerm( f );


    // BC Definition
    //
    UInt dim_fluid = ns->uDof().numTotalDof();
    Vector vec_press(dim_fluid);
    BCVector bcvec(vec_press,dim_fluid,1);

    BCh_u.addBC("Wall",   2, Essential, Full, u_wall,  3);
    BCh_u.addBC("Wall-inflow",   4, Essential, Full, u_wall,  3);
    BCh_u.addBC("Wall-outflow",   5, Essential, Full, u_wall,  3);
    BCh_u.addBC("OutFlow", 3, Natural,  Full, out_flow, 3);


    // Initialization
    //
    std::ofstream outfile;

    Real dt = ns->timestep();
    Real startT = ns->inittime();
    Real T  = ns->endtime();

    outfile.open("res_Q.m", std::ios::app);
    outfile << "xx=[0:" << dt << ":" << T << "]; " << std::endl;
    outfile.close();

    outfile.open("res_DP.m", std::ios::app);
    outfile << "xx=[0:" << dt << ":" << T << "]; " << std::endl;
    outfile.close();


    if(startT > 0.0){
        std::cout << "initialize velocity and pressure with data from file" << std::endl;
        std::ostringstream indexin;
        std::string vinname, cinname;
        indexin << (startT*100);
        vinname = "fluid.res"+indexin.str();
        ns->initialize(vinname);}
    else{
        std::cout << "initialize velocity and pressure with u0 and p0" << std::endl;
        ns->initialize(u0,p0,0.0,dt);
    }

    zeroDModelSolver network(data_file);

    // Temporal loop
    //
    switch(strategy)
    {
        case 1:
        {
            std::cout << "Mean Pressure problem \n";
            std::cout <<" -----------------------------------------------\n";

            BCh_u.addBC("InFlow", 1, Natural,  Full, bcvec , 3);
            Real Q;

            for (Real time=startT+dt ; time <= T; time+=dt)
            {

                Q = ns->flux(1);
                // pressure coming from the 0D Model
                vec_press = ScalarVector( vec_press.size(), -network.getPressureFromQ(time,Q) );

                ns->timeAdvance(f,time);
                ns->iterate(time);

                ns->postProcess();
            }
        }
        break;

        case 2:
        {
            std::cout << "Flow Rate problem \n";
            std::cout <<" -----------------------------------------------\n";

            BCFunctionBase in_flow(uo);
            BCh_u.addBC("InFlow", 1, Natural,   Full, in_flow, 3);

            Real deltaP;

            NavierStokesWithFlux<NS> ns_with_flux(ns);
            ns_with_flux.setFlux( 1, flux_adaptor( 0 ) );
            ns_with_flux.initialize(u0,p0,0.0,dt);

            toEnsight EnsightFilter;
            ns_with_flux.doOnIterationFinish( EnsightFilter  );
            toMedit MeditFilter( data_file );
            ns_with_flux.doOnIterationFinish( MeditFilter  );

            for (Real time=startT+dt ; time <= T; time+=dt){

                deltaP = ns_with_flux.pressure();

                ns_with_flux.setFlux(1, flux_adaptor( network.getQFromPressure(time, deltaP) ) );

                ns_with_flux.iterate( time );
            }
        }
        break;

        default:
            ERROR_MSG("No coupling strategy defined \n");
    }

    outfile.open("res_Q.m", std::ios::app);
    outfile << "    ]; " << std::endl;
    outfile << "figure;" << std::endl;
    outfile << "plot(xx,Q);" << std::endl;
    outfile.close();

    outfile.open("res_DP.m", std::ios::app);
    outfile << "     ]; " << std::endl;
    outfile << "figure;" << std::endl;
    outfile << "plot(xx,DP);" << std::endl;
    outfile.close();


    return 0;
}
