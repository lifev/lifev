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


    //========================================================================================
    //  TEMPORAL LOOP
    //========================================================================================

    operFS oper(data_file);

    UInt maxpf  = 100;
    Real dt     = oper.fluid().timestep();
    Real T      = oper.fluid().endtime();


    Real abstol = 1.e-7;
    Real reltol = 0.;
    Real etamax = 1.e-3;

    int status;
    int maxiter;
    int linesearch = 0;

    std::ofstream nout("num_iter");
    ASSERT(nout,"Error: Output file cannot be opened.");

    UInt dim_solid = oper.solid().dDof().numTotalDof();
    UInt dim_fluid = oper.fluid().uDof().numTotalDof();

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
        oper.fluid().timeAdvance(f,time);
        oper.solid().timeAdvance(f,time);
        oper.setTime(time);

        // displacement prediction

        disp   = 0.;//solid.d();// + dt*(1.5*solid.w() - 0.5*velo_1);

        velo_1 = oper.solid().w();

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

            oper.fluid().postProcess();
            oper.solid().postProcess();
        }
    }

    return 0;
}

