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


    //========================================================================================
    //  TEMPORAL LOOP
    //========================================================================================

    int  method   = data_file("problem/method"  , 0);
    Real defOmega = data_file("problem/defOmega", .01);

    std::cout << std::endl;
    std::cout << "Fluid/Structure interactions";
    std::auto_ptr<operFS> p_oper;
    switch(method)
    {
        case 0:
            std::cout << " -- Fixed Point method" << std::endl;
            std::cout << "--------------------------------------------------"
                      << std::endl;
            std::cout << std::endl;
            p_oper.reset(new fixedPoint(data_file));
            break;
        case 1:
            std::cout << " -- SteklovPoincare method" << std::endl;
            std::cout << "------------------------------------------------------"
                      << std::endl;
            std::cout << std::endl;
            p_oper.reset(new steklovPoincare(data_file));
            break;
        case 2:
            std::cout << " -- exactJacobian method" << std::endl;
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << std::endl;
            p_oper.reset(new exactJacobian(data_file));
            break;
        default:
            p_oper.reset(new steklovPoincare(data_file));
    }

    operFS &oper = *p_oper;

    oper.fluid().showMe();
    oper.solid().showMe();

    UInt maxpf  = 100;
    Real dt     = oper.fluid().timestep();
    Real T      = oper.fluid().endtime();

    Real abstol = 1.e-7;
    Real reltol = 1.e-4;
    Real etamax = 1.e-3;

    int status;
    int maxiter;
    int linesearch = 0;

    std::ofstream nout("num_iter");
    ASSERT(nout,"Error: Output file cannot be opened.");

    UInt dim_solid = oper.solid().dDof().numTotalDof();
    UInt dim_fluid = oper.fluid().uDof().numTotalDof();

    Vector disp(3*dim_solid);
    disp   = ZeroVector( disp.size() );

    Vector velo_1(3*dim_solid);
    velo_1 = ZeroVector( velo_1.size() );

    std::ofstream out_iter("iter");
    std::ofstream out_res ("res");

    oper.setUpBC(fZero, u2);

    oper.fluid().initialize(u0);
    oper.solid().initialize(d0,w0);

    //
    // Temporal loop
    //

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
            status = newton(disp,oper, norm_inf_adaptor(), abstol, reltol, maxiter, etamax,linesearch,out_res,time);
        }
        else
        {
            status = nonLinRichardson(disp, oper, norm_inf_adaptor(), abstol, reltol,
                                      maxiter, etamax, linesearch, out_res,
                                      time, defOmega);
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

    return 0;
}

