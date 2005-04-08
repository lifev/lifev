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
#include <life/lifecore/life.hpp>
#include <life/lifesolver/NavierStokesAleSolverPC.hpp>
#include <life/lifesolver/VenantKirchhofSolver.hpp>
#include <life/lifealg/nonLinRichardson.hpp>
#include <life/lifesolver/fixedPointBase.hpp>
#include <life/lifesolver/exactJacobianBase.hpp>
#include <life/lifefem/dofInterface3Dto3D.hpp>
#include "ud_functions.hpp"
#include <life/lifefem/regionMesh3D_ALE.hpp>
#include "zeroDModelSolver.hpp"


#include <life/lifefilters/ensight7Writer.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>

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

class disp_adaptor
{
public:

    disp_adaptor( Real __displacement)
        :
        _M_displacement( __displacement)
    {}

    Real operator()(Real __time, Real __x, Real __y, Real __z, ID __i)
    {
        switch(__i) {
        case 1:
            return _M_displacement*__x;
            break;
        case 2:
            return _M_displacement*__y;
            break;
        case 3:
            return 0.;
            break;
        default:
            ERROR_MSG("This entrie is not allowed: disp_adatptor");
            break;
        }
    }
private:
    Real  _M_displacement;
};

}

/*

   This programs couples the fluid-structure coupling with a lumped parameter model.
   One strategy is applied:
   Mean pressure : the 3D Model gives the flux information to the 0D Model whereas
                   the 0D Model passes the pressure to the 3D.

*/


int main(int argc, char** argv)
{
#if 0
    using namespace LifeV;


    // Reading from data file
    //
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data_FSI0d", 2, "-f","--file");
    GetPot data_file(data_file_name);

    // Number of boundary conditions for the fluid velocity,
    // solid displacement, and fluid mesh motion
    //
    BCHandler BCh_u;
    BCHandler BCh_d;
    BCHandler BCh_mesh;


    //========================================================================================
    // FLUID AND SOLID SOLVERS
    //========================================================================================
    //
    // The NavierStokes ALE solver
    //
    NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> > fluid(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt,
                                                                   quadRuleTria3pt, quadRuleTetra64pt, quadRuleTria3pt,
                                                                   BCh_u,BCh_mesh);

    // The structural solver
    //
    VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> > solid(data_file, feTetraP1, quadRuleTetra4pt,
                                                                quadRuleTria3pt, BCh_d);

    // The 0d Model
    //
    zeroDModelSolver network(data_file);

    // Outputs
    fluid.showMe();
    solid.showMe();

    UInt dim_solid = solid.dDof().numTotalDof();
    UInt dim_fluid = fluid.uDof().numTotalDof();


    //========================================================================================
    //  DATA INTERFACING BETWEEN BOTH SOLVERS
    //========================================================================================
    //
    // Passing data from the fluid to the structure: fluid load at the interface
    //
    DofInterface3Dto3D dofFluidToStructure(feTetraP1, solid.dDof(), feTetraP1bubble, fluid.uDof());
    dofFluidToStructure.update(solid.mesh(), 1, fluid.mesh(), 1, 0.0);
    BCVectorInterface g_wall(fluid.residual(), dim_fluid, dofFluidToStructure);


    // Passing data from structure to the fluid mesh: motion of the fluid domain
    //
    DofInterface3Dto3D dofStructureToFluidMesh(fluid.mesh().getRefFE(), fluid.dofMesh(),
                                               feTetraP1, solid.dDof());
    dofStructureToFluidMesh.update(fluid.mesh(), 1, solid.mesh(), 1, 0.0);
    BCVectorInterface displ(solid.d(), dim_solid, dofStructureToFluidMesh);

    //========================================================================================
    //  BOUNDARY CONDITIONS
    //========================================================================================
    //

    BCFunctionBase area;
    BCFunctionBase bcf(fZero);

    std::vector<ID> compz(1), compxy(2);
    compz[0]=3;
    compxy[0]=1;
    compxy[1]=2;

    // Boundary conditions for the harmonic extension of the interface solid displacement
    //BCh_mesh.addBC("Interface",            1, Essential, Full,      displ, 3);
    BCh_mesh.addBC("Interface_edge_top",  10, Essential, Component, displ, compxy);
    BCh_mesh.addBC("Interface_edge_top",  10, Essential, Component, bcf,   compz);
    BCh_mesh.addBC("Interface_edge_base", 20, Essential, Component, displ, compxy);
    BCh_mesh.addBC("Interface_edge_base", 20, Essential, Component, bcf,   compz);
    BCh_mesh.addBC("Top",                  3, Essential, Component, bcf,   compz);
    BCh_mesh.addBC("Base",                 2, Essential, Component, bcf,   compz);

    // Boundary conditions for the fluid velocity
    Vector vectorPressure(dim_fluid);
    BCVector bcvecPressure(vectorPressure, dim_fluid, 1);
    BCVector u_wall(fluid.wInterpolated(), dim_fluid);   // Passing w -> u at the interface
    BCh_u.addBC("Wall",       1, Essential, Full, u_wall,        3);
    BCh_u.addBC("Edge_top",  10, Essential, Full, u_wall,        3);
    BCh_u.addBC("Edge_base", 20, Essential, Full, u_wall,        3);
    BCh_u.addBC("InFlow",     2, Natural,   Full, bcvecPressure, 3);


    // Boundary conditions for the solid displacement
    //BCh_d.addBC("Interface", 1, Natural, Full, g_wall, 3);
    BCh_d.addBC("Top",       3, Essential, Component, bcf,compz);
    BCh_d.addBC("Base",      2, Essential, Component, bcf, compz);
    BCh_d.addBC("Base_edge", 4, Essential, Full, area,  3);

    //========================================================================================
    //  COUPLING FS
    //========================================================================================
    //

    UInt method    = data_file("problem/method"    , 0);
    UInt maxpf     = data_file("problem/maxSubIter", 300);
    UInt precond   = data_file("problem/precond"   , 1);
    Real defOmega  = data_file("problem/defOmega"  , 0.01);
    Real abstol    = data_file("problem/abstol"  , 1.e-07);
    Real reltol    = data_file("problem/reltol"  , 1.e-04);
    Real etamax    = data_file("problem/etamax"  , 1.e-03);
    int linesearch = data_file("problem/linesearch"  , 0);

    std::cout << "\n Fluid/Structure interactions \n";
    std::auto_ptr<operFS> __operFSI;
    switch(method)
    {
        case 0:
            std::cout << " -- Fixed Point method \n";
            std::cout << "-------------------------------------------------- \n \n";
            __operFSI.reset(new fixedPoint(fluid,
                                           solid,
                                           data_file,
                                           BCh_u,
                                           BCh_d,
                                           BCh_mesh));
            break;
//         case 1:
//             std::cout << " -- SteklovPoincare method \n";
//             std::cout << "------------------------------------------------------ \n\n";
//             __operFSI.reset(new steklovPoincare(fluid,
//                                                 solid,
//                                                 data_file,
//                                                 BCh_u,
//                                                 BCh_d,
//                                                 BCh_mesh));
//             if (precond == 2) defOmega = -1.;
//             break;
        case 2:
            std::cout << " -- exactJacobian method\n";
            std::cout << "----------------------------------------------------\n\n";
            __operFSI.reset(new exactJacobian(fluid,
                                              solid,
                                              data_file,
                                              BCh_u,
                                              BCh_d,
                                              BCh_mesh));
            break;
        default:
            ERROR_MSG("No FSI coupling strategy defined \n");
    }

    operFS &operFSI = *__operFSI;


    //========================================================================================
    //  TEMPORAL LOOP
    //========================================================================================
    //

    int status;
    int maxiter;

    Real dt = fluid.timestep();
    Real T  = fluid.endtime();

    fluid.initialize(u0);
    solid.initialize(d0,w0);

    std::ofstream nout("num_iter");
    ASSERT(nout,"Error: Output file cannot be opened.");

    Vector disp(3*dim_solid);
    disp = ZeroVector( disp.size() );

    Vector velo_1(3*dim_solid);
    velo_1 = ZeroVector( velo_1.size() );

    std::ofstream out_iter("iter");
    std::ofstream out_res("res");

    Real Qin;
    Real Qout;
    Real deltaP;
    Real displacement;

    std::ofstream _M_outfile;

    _M_outfile.open("res_Qin.m", std::ios::app);
    _M_outfile << "xx=[0:" << dt << ":" << T << "]; " << std::endl;
    _M_outfile << "Qin = [ " << std::endl;
    _M_outfile.close();

    _M_outfile.open("res_Qout.m", std::ios::app);
    _M_outfile << "xx=[0:" << dt << ":" << T << "]; " << std::endl;
    _M_outfile << "Qout = [ " << std::endl;
    _M_outfile.close();

    _M_outfile.open("res_DP.m", std::ios::app);
    _M_outfile << "xx=[0:" << dt << ":" << T << "]; " << std::endl;
    _M_outfile << "DP = [ " << std::endl;
    _M_outfile.close();

    // Temporal loop
    //

    Chrono chrono;
    chrono.start();

    for (Real time=dt; time <= T; time+=dt) {

        Qin    = fluid.flux(2);
        Qout   = fluid.flux(3);
        deltaP = network.getPressureFromQFSI(time, Qin, Qout);

        vectorPressure = ScalarVector(dim_fluid, -deltaP);

        displacement = network.radius()-1;

        area.setFunction(disp_adaptor(displacement));
        BCh_d.modifyBC( 4, area);

        _M_outfile.open("res_Qin.m", std::ios::app);
        _M_outfile << "     " << -Qin << ";" << std::endl;
        _M_outfile.close();

        _M_outfile.open("res_Qout.m", std::ios::app);
        _M_outfile << "     " << Qout << ";" << std::endl;
        _M_outfile.close();

        _M_outfile.open("res_DP.m", std::ios::app);
        _M_outfile << "      " <<  deltaP << ";" << std::endl;
        _M_outfile.close();

        // displacement prediction
        disp = solid.d() + dt*(1.5*solid.w() - 0.5*velo_1);

        velo_1 = solid.w();

        std::cout << "norm( disp ) init = " << norm_inf(disp) << std::endl;
        std::cout << "norm( velo_! ) init = " << norm_inf(velo_1) << std::endl;

        maxiter = maxpf;

        // the newton solver
        if (method == 2)
        {
            status = newton(disp, operFSI, norm_inf_adaptor(),
                            abstol, reltol, maxiter, etamax,
                            linesearch, out_res, time);
        }
        else
        {
            status = nonLinRichardson(disp, operFSI, norm_inf_adaptor(),
                                      abstol, reltol, maxiter, etamax,
                                      linesearch, out_res, time, defOmega);
        }


        if(status == 1) {
            std::cout << "Inners iterations failed\n";
            exit(1);
        }
        else {
            std::cout << "End of time "<< time << std::endl;
            std::cout << "Number of inner iterations       : " << maxiter << std::endl;
            out_iter << time << " " << maxiter << " " << operFSI.nbEval() << std::endl;
            fluid.postProcess();
            solid.postProcess();
        }

    }

    chrono.stop();
    out_res << "Total computation time = " << chrono.diff() << "\n";

    _M_outfile.open("res_Qin.m", std::ios::app);
    _M_outfile << "    ]; " << std::endl;
    _M_outfile << "figure;" << std::endl;
    _M_outfile << "plot(xx,Qin);" << std::endl;
    _M_outfile.close();

    _M_outfile.open("res_Qout.m", std::ios::app);
    _M_outfile << "    ]; " << std::endl;
    _M_outfile << "figure;" << std::endl;
    _M_outfile << "plot(xx,Qout);" << std::endl;
    _M_outfile.close();

    _M_outfile.open("res_DP.m", std::ios::app);
    _M_outfile << "     ]; " << std::endl;
    _M_outfile << "figure;" << std::endl;
    _M_outfile << "plot(xx,DP);" << std::endl;
    _M_outfile.close();
#endif
    return 0;
}

