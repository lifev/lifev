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
#include "steklovPoincareBase.hpp"
#include "fixedPointBase.hpp"
#include "exactJacobianBase.hpp"
#include "dofInterface3Dto3D.hpp"
#include "ud_functions.hpp"
#include "regionMesh3D_ALE.hpp"


/*

This programs couples the Navier-Stokes and (linear) Elastodynamic equations
At each time step the resulting non-linear coupled problem is solved via
Domain Decomposition method.

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

    // Outputs
    fluid.showMe();
    solid.showMe();

    UInt method   = data_file("problem/method"    , 0);
    UInt maxpf    = data_file("problem/maxSubIter", 300);
    UInt precond  = data_file("problem/precond"   , 1);
    Real defOmega = data_file("problem/defOmega"  , 0.01);


    std::cout << std::endl;
    std::cout << "Fluid/Structure interactions";
    std::cout << std::scientific;
    std::auto_ptr<operFS> p_oper;
    switch(method)
    {
        case 0:
            std::cout << " -- Fixed Point method" << std::endl;
            std::cout << "--------------------------------------------------"
                      << std::endl;
            std::cout << std::endl;
            p_oper.reset(new fixedPoint(fluid,
                                        solid,
                                        data_file,
                                        BCh_u,
                                        BCh_d,
                                        BCh_mesh));
            break;
        case 1:
            std::cout << " -- SteklovPoincare method" << std::endl;
            std::cout << "------------------------------------------------------"
                      << std::endl;
            std::cout << std::endl;
            p_oper.reset(new steklovPoincare(fluid,
                                             solid,
                                             data_file,
                                             BCh_u,
                                             BCh_d,
                                             BCh_mesh));
            if (precond == 2) defOmega = -1.;
            break;
        case 2:
            std::cout << " -- exactJacobian method" << std::endl;
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << std::endl;
            p_oper.reset(new exactJacobian(fluid,
                                           solid,
                                           data_file,
                                           BCh_u,
                                           BCh_d,
                                           BCh_mesh));
            break;
        default:
            p_oper.reset(new steklovPoincare(fluid,
                                             solid,
                                             data_file,
                                             BCh_u,
                                             BCh_d,
                                             BCh_mesh));
    }

    operFS &oper = *p_oper;

    UInt dim_solid = oper.solid().dDof().numTotalDof();
    UInt dim_fluid = oper.fluid().uDof().numTotalDof();


    //========================================================================================
    //  DATA INTERFACING BETWEEN BOTH SOLVERS
    //========================================================================================
    //
    // Passing data from the fluid to the structure: fluid load at the interface
    //
    DofInterface3Dto3D dofFluidToStructure(feTetraP1, oper.solid().dDof(), feTetraP1bubble, oper.fluid().uDof());
    dofFluidToStructure.update(oper.solid().mesh(), 1, oper.fluid().mesh(), 1, 0.0);
    BCVectorInterface g_wall(oper.fluid().residual(), dim_fluid, dofFluidToStructure);


    // Passing data from structure to the fluid mesh: motion of the fluid domain
    //
    DofInterface3Dto3D dofStructureToFluidMesh(oper.fluid().mesh().getRefFE(), oper.fluid().dofMesh(),
                                               feTetraP1, oper.solid().dDof());
    dofStructureToFluidMesh.update(oper.fluid().mesh(), 1, oper.solid().mesh(), 1, 0.0);
    BCVectorInterface displ(oper.solid().d(), dim_solid, dofStructureToFluidMesh);



    // Passing data from structure to the fluid: solid velocity at the interface velocity
    //
    DofInterface3Dto3D dofMeshToFluid(feTetraP1bubble, oper.fluid().uDof(), feTetraP1bubble, oper.fluid().uDof() );
    dofMeshToFluid.update(oper.fluid().mesh(), 1, oper.fluid().mesh(), 1, 0.0);
    BCVectorInterface u_wall(oper.fluid().wInterpolated(), dim_fluid,dofMeshToFluid);

    //========================================================================================
    //  BOUNDARY CONDITIONS
    //========================================================================================
    //
    // Boundary conditions for the harmonic extension of the
    // interface solid displacement
    BCFunctionBase bcf(fZero);
//    BCh_mesh.addBC("Interface", 1, Essential, Full, displ, 3);
    BCh_mesh.addBC("Top",       3, Essential, Full, bcf,   3);
    BCh_mesh.addBC("Base",      2, Essential, Full, bcf,   3);
    BCh_mesh.addBC("Edges",    20, Essential, Full, bcf,   3);

    // Boundary conditions for the fluid velocity
    BCFunctionBase in_flow(u2);
    BCh_u.addBC("Wall",   1,  Essential, Full, u_wall,  3);
    BCh_u.addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
    BCh_u.addBC("Edges",  20, Essential, Full, bcf,     3);

    // Boundary conditions for the solid displacement
//    BCh_d.addBC("Interface", 1, Natural, Full, g_wall, 3);
    BCh_d.addBC("Top",       3, Essential, Full, bcf,  3);
    BCh_d.addBC("Base",      2, Essential, Full, bcf,  3);


    //========================================================================================
    //  TEMPORAL LOOP
    //========================================================================================



    oper.fluid().showMe();
    oper.solid().showMe();

    Real dt     = oper.fluid().timestep();
    Real T      = oper.fluid().endtime();

//     Real abstol = 1.e-7;
//     Real reltol = 1.e-4;
//     Real etamax = 1.e-3;

    Real abstol = data_file("problem/abstol"  , 1.e-07);
    Real reltol = data_file("problem/reltol"  , 1.e-04);
    Real etamax = data_file("problem/etamax"  , 1.e-03);

    int status;
    int maxiter;

    int linesearch = data_file("problem/linesearch"  , 0);

    std::ofstream nout("num_iter");
    ASSERT(nout,"Error: Output file cannot be opened.");

    Vector disp(3*dim_solid);
    disp   = ZeroVector( disp.size() );

    Vector velo_1(3*dim_solid);
    velo_1 = ZeroVector( velo_1.size() );

    std::ofstream out_iter("iter");
    std::ofstream out_res ("res");

    oper.fluid().initialize(u0);
    oper.solid().initialize(d0,w0);

    //
    // Temporal loop
    //

    Chrono chrono;
    chrono.start();

    for (Real time=dt; time <= T; time+=dt)
    {
        oper.fluid().timeAdvance(f,time);
        oper.solid().timeAdvance(f,time);
        oper.setTime(time);

        // displacement prediction

        disp   = oper.solid().d() + dt*(1.5*oper.solid().w() - 0.5*velo_1);

        velo_1 = oper.solid().w();

        std::cout << "norm( disp   ) init = " << norm_inf(disp)   << std::endl;
        std::cout << "norm( velo_1 ) init = " << norm_inf(velo_1) << std::endl;

        maxiter = maxpf;

        // the newton solver

        if (method == 2)
        {
            status = newton(disp, oper, norm_inf_adaptor(),
                        abstol, reltol, maxiter, etamax,
                        linesearch, out_res, time);
        }
        else
        {
            status = nonLinRichardson(disp, oper, norm_inf_adaptor(),
                                      abstol, reltol, maxiter, etamax,
                                      linesearch, out_res, time, defOmega);
        }

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

            oper.fluid().postProcess();
            oper.solid().postProcess();
        }
    }

    chrono.stop();
    out_res << "Total computation time = " << chrono.diff()
            << std::endl;
    return 0;
}

