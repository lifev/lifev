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
#include "life.hpp"
#include "NavierStokesAleSolverPC.hpp"
#include "VenantKirchhofSolver.hpp"
#include "operFS.hpp"
#include "dofInterface3Dto3D.hpp"
#include "ud_functions.hpp"
#include "regionMesh3D_ALE.hpp"


/*

   This programs couples the Navier-Stokes and (linear) Elastodynamic equations
   At each time step the resulting non-linear coupled problem is solved using a
   quasi Newton method. The fluid jacobian is approximated using a reduced model
   (see Gerbeau-Vidrascu M2AN 2003)

   The present test simulates the pressure wave propagation in a curved cylindrical vessel

*/


int main(int argc, char** argv)
{
    using namespace LifeV;
    using namespace std;
    typedef  boost::shared_ptr<DofInterface3Dto3D> dof_interface_type;


    // Reading from data file
    //
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);




    // solid displacement, fluid velocity and mesh displacement BC's
    //
    BCHandler BCh_u;
    BCHandler BCh_d;
    BCHandler BCh_mesh;


    // fluid solver
    //
    NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> > fluid(data_file,
								   feTetraP1bubble,
								   feTetraP1,
								   quadRuleTetra64pt,
                                                                   quadRuleTria3pt,
								   quadRuleTetra64pt,
								   quadRuleTria3pt,
                                                                   BCh_u,BCh_mesh);

    // structural solver
    //
    VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> > solid(data_file,
								feTetraP1,
								quadRuleTetra4pt,
                                                                quadRuleTria3pt,
								BCh_d);


    // Outputs
    fluid.showMe();
    solid.showMe();

    UInt dim_solid = solid.dDof().numTotalDof();
    UInt dim_fluid = fluid.uDof().numTotalDof();
    UInt dim_reducedfluid = fluid.pDof().numTotalDof();


    // passing data from the fluid to the structure: fluid load at the interface
    //
    dof_interface_type  dofFluidToStructure( new DofInterface3Dto3D(feTetraP1,
								    solid.dDof(),
								    feTetraP1bubble,
								    fluid.uDof()) );
    dofFluidToStructure->update(solid.mesh(), 1, fluid.mesh(), 1, 0.0);

    // passing data from structure to the fluid mesh: motion of the fluid domain
    //
    dof_interface_type dofStructureToFluidMesh( new  DofInterface3Dto3D(fluid.mesh().getRefFE(),
									fluid.dofMesh(),
									feTetraP1,
									solid.dDof()) );
    dofStructureToFluidMesh->update(fluid.mesh(), 1, solid.mesh(), 1, 0.0);


    // BC's for the harmonic extension of the
    // interface solid displacement
    BCFunctionBase bcf(fZero);
    BCVectorInterface displ(solid.d(), dim_solid, dofStructureToFluidMesh );
    BCh_mesh.addBC("Interface", 1, Essential, Full, displ, 3);
    BCh_mesh.addBC("Top",       3, Essential, Full, bcf,   3);
    BCh_mesh.addBC("Base",      2, Essential, Full, bcf,   3);
    BCh_mesh.addBC("Edges",    20, Essential, Full, bcf,   3);

    // BC's for the fluid velocity u
    BCFunctionBase in_flow(u2);
    BCVector u_wall(fluid.wInterpolated(), dim_fluid);   // Passing w -> u at the interface
    BCh_u.addBC("Wall",        1, Essential, Full, u_wall,  3);
    BCh_u.addBC("Wall_Edges", 20, Essential, Full, u_wall,  3);
    BCh_u.addBC("InFlow",      2, Natural,   Full, in_flow, 3);

    // BC's for the solid displacement d
    BCVectorInterface g_wall(fluid.residual(), dim_fluid,  dofFluidToStructure );
    BCh_d.addBC("Interface", 1, Natural, Full, g_wall, 3);
    BCh_d.addBC("Top",       3, Essential, Full, bcf,  3);
    BCh_d.addBC("Base",      2, Essential, Full, bcf,  3);



    // pressure and displacement variations BC's
    BCHandler BCh_dp; // approximated using a reduced model (see Gerbeau-Vidrascu M2AN 2003)
    BCHandler BCh_dz;

    // the coupled FSI operator
    operFS oper(fluid,solid,BCh_dp,BCh_dz,data_file);

    // passing data bettwen the reduced linearized fluid to the linearized solver solid:
    // reduced fluid pressure
    dof_interface_type dofReducedFluidToStructure( new DofInterface3Dto3D(feTetraP1,
									  solid.dDof(),
									  feTetraP1,
									  fluid.pDof()) );
    dofReducedFluidToStructure->update(solid.mesh(), 1, fluid.mesh(), 1, 0.0);

    // solid acceleration
    dof_interface_type dofStructureToReducedFluid( new DofInterface3Dto3D(feTetraP1,
									  fluid.pDof(),
									  feTetraP1,
									  solid.dDof()) );
    dofStructureToReducedFluid->update(fluid.mesh(), 1, solid.mesh(), 1, 0.0);


    // Boundary conditions for dp
    BCVectorInterface da_wall(oper.da(), dim_solid, dofStructureToReducedFluid,2); // type  = 2
    BCh_dp.addBC("Wall",        1, Natural,   Scalar, da_wall);
    BCh_dp.addBC("Wall_Edges", 20, Essential, Scalar, bcf);
    BCh_dp.addBC("InFlow",      2, Essential, Scalar, bcf);
    BCh_dp.addBC("OutFlow",     3, Essential, Scalar, bcf);


    // Boundary conditions for dz
    BCVectorInterface dg_wall(oper.minusdp(), dim_reducedfluid, dofReducedFluidToStructure, 1); // type = 1
    BCh_dz.addBC("Interface", 1, Natural,   Full, dg_wall, 3);
    BCh_dz.addBC("Top",       3, Essential, Full, bcf,     3);
    BCh_dz.addBC("Base",      2, Essential, Full, bcf,     3);



    //  TEMPORAL LOOP

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
        disp = solid.d() + dt*(1.5*solid.w() - 0.5*velo_1);

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
            fluid.postProcess();
            solid.postProcess();
        }
    }
    return 0;
}
