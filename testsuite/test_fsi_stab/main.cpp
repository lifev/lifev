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
#include <life/lifesolver/NavierStokesAleSolver.hpp>
#include <life/lifesolver/VenantKirchhofSolver.hpp>
#include "operFS.hpp"
#include <life/lifefem/dofInterface3Dto3D.hpp>
#include "ud_functions.hpp"
#include <life/lifefem/regionMesh3D_ALE.hpp>
#include <life/lifefilters/medit.hpp>

/*

   This programs couples the Navier-Stokes and (linear) Elastodynamic equations
   At each time step the resulting non-linear coupled problem is solved via
   Newton's method. The jacobian is fully computed (see Fernandez-Moubachir 2003,2004 )

   The present test simulates the pressure wave propagation in a curved cylindrical vessel

*/


int main(int argc, char** argv)
{
    using namespace LifeV;
    using namespace std;



    // Reading from data file
    //
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);



    // Number of boundary conditions for the fluid velocity,
    // solid displacement, and fluid mesh motion
    //
    BCHandler bcHu; 
    BCHandler bcHd; 
    BCHandler bcHmesh; 


    //========================================================================================
    // FLUID AND SOLID SOLVERS
    //========================================================================================
    //
    // The NavierStokes ALE solver
    //
    NavierStokesAleSolver< RegionMesh3D_ALE<LinearTetra> >   fluid(data_file, feTetraP1, quadRuleTetra4pt,
								     quadRuleTria3pt, bcHu,bcHmesh);

    // The structural solver
    //
    VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> > solid(data_file, feTetraP1, quadRuleTetra4pt,
                                                                quadRuleTria3pt, bcHd);

    Medit medit(fluid);

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
    boost::shared_ptr<DofInterface3Dto3D> dofFluidToStructure( new DofInterface3Dto3D(feTetraP1, solid.dDof(), feTetraP1, fluid.uDof()) );
    dofFluidToStructure->update(solid.mesh(), 1, fluid.mesh(), 1, 0.0);
    BCVectorInterface g_wall(fluid.residual(), dim_fluid,  dofFluidToStructure );


    // Passing data from structure to the fluid mesh: motion of the fluid domain
    //
    boost::shared_ptr<DofInterface3Dto3D> dofStructureToFluidMesh( new  DofInterface3Dto3D(fluid.mesh().getRefFE(), fluid.uDof(),
											    feTetraP1, solid.dDof()) );
    dofStructureToFluidMesh->update(fluid.mesh(), 1, solid.mesh(), 1, 0.0);
    BCVectorInterface displ(solid.disp(), dim_solid, dofStructureToFluidMesh );


    // Passing data from structure to the fluid pressure: reduced fluid
    //
    boost::shared_ptr<DofInterface3Dto3D> dofStructureToFluid( new  DofInterface3Dto3D(feTetraP1, fluid.pDof(),
										       feTetraP1, solid.dDof()) );
    dofStructureToFluid->update(fluid.mesh(), 1, solid.mesh(), 1, 0.0);
  
 

    //========================================================================================
    //  BOUNDARY CONDITIONS
    //========================================================================================
    //
    // Boundary conditions for the harmonic extension of the
    // interface solid displacement
    BCFunctionBase bcf(fZero);
    bcHmesh.addBC("Interface", 1, Essential, Full, displ, 3);
    bcHmesh.addBC("Top",       3, Essential, Full, bcf,   3);
    bcHmesh.addBC("Base",      2, Essential, Full, bcf,   3);
    bcHmesh.addBC("Edges",    20, Essential, Full, bcf,   3);

    // Boundary conditions for the fluid velocity
    BCFunctionBase in_flow(u2);  
    BCVector u_wall(fluid.w(), dim_fluid);   // Passing w -> u at the interface
    bcHu.addBC("Wall",        1, Essential, Full, u_wall,  3);    
    bcHu.addBC("Wall_Edges", 20, Essential, Full, u_wall,  3);
    bcHu.addBC("InFlow",      2, Natural,   Full, in_flow, 3);
    
    // Boundary conditions for the solid displacement
    bcHd.addBC("Interface", 1, Natural, Full, g_wall, 3);
    bcHd.addBC("Top",       3, Essential, Full, bcf,  3);
    bcHd.addBC("Base",      2, Essential, Full, bcf,  3);


    //========================================================================================
    //  COUPLED FSI OPERATOR
    //========================================================================================
    //
    //
    BCHandler bcHdu; 
    BCHandler bcHdp;
    BCHandler bcHdz; 
    BCHandler sum_dp;

    operFS<  NavierStokesAleSolver< RegionMesh3D_ALE<LinearTetra> > ,  
    VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >  > 
      oper(fluid, solid, bcHdu, bcHdp, sum_dp, bcHdz, data_file);

    // Passing data from fluid to the structure: du -> dz
    //
    BCVectorInterface dg_wall(fluid.residual(), dim_fluid, dofFluidToStructure );
    BCVector dp_wall( fluid.dp(),          dim_fluid, 1);

    // Boundary conditions for du
    BCVector du_wall(fluid.dw(), dim_fluid); // dw -> du
    bcHdu.addBC("Wall",         1, Essential, Full, du_wall, 3);
    bcHdu.addBC("Wall_Edges",  20, Essential, Full, du_wall, 3);
    bcHdu.addBC("Inlet", 2, Natural, Full,  bcf, 3);

    // Boundary conditions for dp
    BCVectorInterface da_wall(oper.da(), dim_solid, dofStructureToFluid, 2); // type  = 2
    bcHdp.addBC("Wall",        1, Natural,   Scalar, da_wall);
    bcHdp.addBC("Wall_Edges", 20, Essential, Scalar, bcf);
    bcHdp.addBC("InFlow",      2, Essential, Scalar, bcf);
    bcHdp.addBC("OutFlow",     3, Essential, Scalar, bcf);

    sum_dp.addBC("Interface", 1, Natural, Full, dp_wall, 3); 

    // Boundary conditions for dz
    bcHdz.addBC("Interface", 1, Natural,   Full, dg_wall, 3);
    bcHdz.addBC("Top",       3, Essential, Full, bcf,     3);
    bcHdz.addBC("Base",      2, Essential, Full, bcf,     3);



    //========================================================================================
    //  TEMPORAL LOOP
    //========================================================================================
    //

    UInt maxpf = 100;
    Real dt = fluid.timestep();
    Real T  = fluid.endtime();
    fluid.initialize(u0);
    solid.initialize(d0,w0);
    Real abstol=1.0e-7, reltol=0.0, etamax=1.e-3;
    int status,maxiter,linesearch=0;

    ofstream nout("num_iter");
    ASSERT(nout,"Error: Output file cannot be opened.");

    Vector disp(3*dim_solid);
    disp = ZeroVector( disp.size() );

    Vector velo_1(3*dim_solid);
    velo_1 = ZeroVector( velo_1.size() );

    ofstream out_iter("iter");
    ofstream out_res("res");

    // Temporal loop
    //
    for (Real time=dt; time <= T; time+=dt) {

        fluid.timeAdvance(f,time);
        solid.timeAdvance(f,time);
        oper.setTime(time);

        // displacement prediction
        disp = solid.disp()+ dt*(1.5*solid.w() - 0.5*velo_1);

        velo_1 = solid.w();

        cout << "norm( disp ) init = " << norm_inf(disp) << endl;
        cout << "norm( velo_! ) init = " << norm_inf(velo_1) << endl;

        maxiter = maxpf;

        // the newton solver
        status = newton(disp,oper, norm_inf_adaptor(), abstol, reltol, maxiter, etamax,linesearch,out_res,time);

        if(status == 1) {
            cout << "Inners iterations failed\n";
            exit(1);
        }
        else {
            cout << "End of time "<< time << endl;
            cout << "Number of inner iterations       : " << maxiter << endl;
            out_iter << time << " " << maxiter << " " << oper.nbEval() << endl;
            medit.postProcessALE(fluid);
            solid.postProcess();
        }
    }
    return 0;
}

