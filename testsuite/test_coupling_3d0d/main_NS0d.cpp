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
#include "chrono.hpp"
#include "ud_functions.hpp"
#include "GetPot.hpp"

#include "meanPressure.hpp"
#include "flowRate.hpp"

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

}
/**
   This test couples the Navier-Stokes equations (3D Model) with a lumped parameter model
   (0D Model - electric network).
   Two couple strategies ares used:

   Mean pressure : the 3D Model gives the flux information to the 0D Model whereas
                   the 0D Model passes the pressure to the 3D.

   Flow rate     : the 3D Model gives the pressure information to the 0D Model whereas
                   the 0D Model passes the flux to the 3D.

*/

int main(int argc, char** argv)
{
    using namespace LifeV;

    // Reading from data file
    //
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data_NS0d", 2, "-f","--file");
    GetPot data_file(data_file_name);

    int strategy = data_file("problem/strategy", 0 );

    // BC Definition
    //
    BCFunctionBase u_wall(u1);
    BCFunctionBase out_flow(u1);
    BCFunctionBase in_flow(uo);
    BCHandler BCh_u(5);

    BCh_u.addBC("Wall",   2, Essential, Full, u_wall,  3);
    BCh_u.addBC("Wall-inflow",   4, Essential, Full, u_wall,  3);
    BCh_u.addBC("Wall-outflow",   5, Essential, Full, u_wall,  3);
    BCh_u.addBC("OutFlow", 3, Natural,  Full, out_flow, 3);
    BCh_u.addBC("InFlow", 1, Natural,   Full, in_flow, 3);

    std::cout << "\n Coupling 3D and 0D models \n ";
    std::auto_ptr<C3d0d> __coupling;

    switch(strategy)
        {
        case 1:
            {
                std::cout << "Mean Pressure problem \n";
                std::cout <<" -----------------------------------------------\n";
                __coupling.reset(new meanPressure(data_file, BCh_u, 1));
            }
            break;

        case 2:
            {
                std::cout << "Flow Rate problem \n";
                std::cout <<" -----------------------------------------------\n";
                __coupling.reset(new flowRate(data_file, BCh_u, 1));

                toEnsight EnsightFilter;
                __coupling->nsFlux().doOnIterationFinish( EnsightFilter  );
                toMedit MeditFilter( data_file );
                __coupling->nsFlux().doOnIterationFinish( MeditFilter  );

            }
            break;

        default:
            ERROR_MSG("No coupling strategy defined \n");
        }

    C3d0d &coupling = *__coupling;

    coupling.ns()->showMe();
    coupling.startWrtMatlab();

    // Initialization
    //
    Real dt     = coupling.ns()->timestep();
    Real startT = coupling.ns()->inittime();
    Real T      = coupling.ns()->endtime();

    std::cout << "initialize velocity and pressure with u0 and p0 \n";
    coupling.initialize(u0,p0,0.0,dt);

    // Temporal loop
    //
    for (Real time=startT+dt ; time <= T; time+=dt)
        {

            coupling.timeAdvance(f,time);
            coupling.iterate(time);

            if (strategy==1)
                coupling.ns()->postProcess();// The observers are not yet implemented in NSSolverPC
                                             // so it is used the "traditional" post proc in this case
        }

    coupling.finishWrtMatlab();


    return 0;
}
