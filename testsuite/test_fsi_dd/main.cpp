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
#include "NavierStokesAleSolverPC.hpp"
#include "VenantKirchhofSolver.hpp"
#include "nonLinRichardson.hpp"
#include "operFS.hpp"
#include "vectorNorms.hpp"
#include "dofInterface3Dto3D.hpp"
#include "ud_functions.hpp"
#include "regionMesh3D_ALE.hpp"


/*

This programs couples the Navier-Stokes and (linear) Elastodynamic equations
At each time step the resulting non-linear coupled problem is solved via
Domain Decomposition method.

The present test simulates the pressure wave propagation in a curved cylindrical vessel

Based on Miguel Fernandez's test_fsi_newton. See:
A Newton method using exact jacobians for solving fluid-structure coupling
Miguel Fernandez, Marwan Moubachir
*/


int main(int argc, char** argv)
{
    using namespace LifeV;

    // Reading from data file
    //
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);



    // Number of boundary conditions for the fluid velocity,
    // solid displacement, and fluid mesh motion

    BC_Handler BCh_u(3);
    BC_Handler BCh_d(3);
    BC_Handler BCh_mesh(4);


    //========================================================================================
    // FLUID AND SOLID SOLVERS
    //========================================================================================

    //
    // The NavierStokes ALE solver
    //
    NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >
        fluid(data_file,
              feTetraP1bubble,
              feTetraP1,
              quadRuleTetra64pt,
              quadRuleTria3pt,
              quadRuleTetra64pt,
              quadRuleTria3pt,
              BCh_u, BCh_mesh);
    //
    // The structural solver
    //
    VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >
        solid(data_file,
              feTetraP1,
              quadRuleTetra4pt,
              quadRuleTria3pt,
              BCh_d);
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

    DofInterface3Dto3D dofFluidToStructure(feTetraP1,
                                           solid.dDof(),
                                           feTetraP1bubble,
                                           fluid.uDof());
    dofFluidToStructure.update(solid.mesh(),
                               1,
                               fluid.mesh(),
                               1,
                               0.);

    BCVector_Interface g_wall(fluid.residual(),
                              dim_fluid,
                              dofFluidToStructure);
    
    //
    // Passing data from structure to the solid mesh: motion of the solid domain
    //
    DofInterface3Dto3D dofStructureToSolid(feTetraP1,
                                           solid.dDof(),
                                           feTetraP1,
                                           solid.dDof() );
    dofStructureToSolid.update(solid.mesh(),
                               1,
                               solid.mesh(),
                               1,
                               0.);

    BCVector_Interface d_wall(solid.d(),
                              dim_solid,
                              dofStructureToSolid);
    //
    // Passing data from structure to the fluid mesh: motion of the fluid domain
    //
    DofInterface3Dto3D dofStructureToFluidMesh(fluid.mesh().getRefFE(),
                                               fluid.dofMesh(),
                                               feTetraP1,
                                               solid.dDof());
    dofStructureToFluidMesh.update(fluid.mesh(),
                                   1,
                                   solid.mesh(),
                                   1,
                                   0.0);
    BCVector_Interface displ(solid.d(),
                             dim_solid,
                             dofStructureToFluidMesh);
    //
    // Passing data from structure to the fluid: solid velocity at the interface velocity
    //
    DofInterface3Dto3D dofMeshToFluid(feTetraP1bubble,
                                      fluid.uDof(),
                                      feTetraP1bubble,
                                      fluid.uDof() );
    dofMeshToFluid.update(fluid.mesh(),
                          1,
                          fluid.mesh(),
                          1,
                          0.0);
    BCVector_Interface u_wall(fluid.wInterpolated(),
                              dim_fluid,
                              dofMeshToFluid);
    //========================================================================================
    //  BOUNDARY CONDITIONS
    //========================================================================================

    // Boundary conditions for the harmonic extension of the
    // interface solid displacement
    BCFunction_Base bcf(fZero);
    BCh_mesh.addBC("Interface", 1, Essential, Full, displ, 3);
    BCh_mesh.addBC("Top",       3, Essential, Full, bcf,   3);
    BCh_mesh.addBC("Base",      2, Essential, Full, bcf,   3);
    BCh_mesh.addBC("Edges",    20, Essential, Full, bcf,   3);

    // Boundary conditions for the fluid velocity
    BCFunction_Base in_flow(u2);
    BCh_u.addBC("Wall",   1,  Essential, Full, u_wall,  3);
    BCh_u.addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
    BCh_u.addBC("Edges",  20, Essential, Full, bcf,     3);

    // Boundary conditions for the solid displacement
    BCh_d.addBC("Interface", 1, Essential, Full, d_wall, 3);
    BCh_d.addBC("Top",       3, Essential, Full, bcf,  3);
    BCh_d.addBC("Base",      2, Essential, Full, bcf,  3);


    //========================================================================================
    //  COUPLED FSI LINEARIZED OPERATORS
    //========================================================================================

    BC_Handler BCh_du(2);
    BC_Handler BCh_dz(3);

    operFS oper(fluid, solid, BCh_du, BCh_dz);

    // Passing the residue to the linearized fluid: \sigma -> du
    //
    // rem: for now: no fluid.dwInterpolated().
    //      In the future this could be relevant

    BCVector_Interface du_wall(oper.residualFSI(),
                               dim_fluid,
                               dofStructureToFluidMesh);
    // Passing the residual to the linearized structure: \sigma -> dz
    BCVector_Interface dg_wall(oper.residualFSI(),
                               dim_fluid,
                               dofFluidToStructure);
    // Boundary conditions for du

    BCh_du.addBC("Wall",   1,  Natural  , Full, du_wall,  3);
    BCh_du.addBC("Edges",  20, Essential, Full, bcf,      3);

    // Boundary conditions for dz

    BCh_dz.addBC("Interface", 1, Natural  , Full, dg_wall, 3);
    BCh_dz.addBC("Top",       3, Essential, Full, bcf,     3);
    BCh_dz.addBC("Base",      2, Essential, Full, bcf,     3);


    //========================================================================================
    //  TEMPORAL LOOP
    //========================================================================================


    UInt maxpf  = 100;
    Real dt     = fluid.timestep();
    Real T      = fluid.endtime();

    fluid.initialize(u0);
    solid.initialize(d0,w0);

    Real abstol = 1.e-6;
    Real reltol = 0.0;
    Real etamax = 1.e-3;

    int status;
    int maxiter;
    int linesearch = 0;

    std::ofstream nout("num_iter");
    ASSERT(nout,"Error: Output file cannot be opened.");

    Vector disp(3*dim_solid);
    disp   = 0.0;

    Vector velo_1(3*dim_solid);
    velo_1 = 0.0;

    std::ofstream out_iter("iter");
    std::ofstream out_res ("res");

    //
    // Temporal loop
    //

    for (Real time=dt; time <= T; time+=dt)
    {
        fluid.timeAdvance(f,time);
        solid.timeAdvance(f,time);
        oper.setTime(time);

        // displacement prediction

        disp   = solid.d() + dt*(1.5*solid.w() - 0.5*velo_1);

        velo_1 = solid.w();

        std::cout << "norm( disp   ) init = " << maxnorm(disp)   << std::endl;
        std::cout << "norm( velo_1 ) init = " << maxnorm(velo_1) << std::endl;

        maxiter = maxpf;

        // the newton solver

        status = nonLinRichardson(disp, oper, maxnorm, abstol, reltol,
                        maxiter, etamax, linesearch, out_res,
                        time, 0.1);
//        status = newton(disp,oper, maxnorm, abstol, reltol, maxiter, etamax,linesearch,out_res,time);

        if(status == 1)
        {
            std::cout << "Inners iterations failed\n";
            exit(1);
        }
        else
        {
            std::cout << "End of time "<< time << std::endl;
            std::cout << "Number of inner iterations       : "
                      << maxiter << std::endl;
            out_iter << time << " " << maxiter << " "
                     << oper.nbEval() << std::endl;

            fluid.postProcess();
            solid.postProcess();
        }
    }

    return 0;
}

