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
#include <life/lifesolver/NavierStokesAleSolverCT.hpp>
#include <life/lifesolver/VenantKirchhofSolver.hpp>
#include "operFS.hpp"
#include <life/lifefem/dofInterface3Dto3D.hpp>
#include "ud_functions.hpp"
#include <life/lifefem/regionMesh3D_ALE.hpp>
#include <life/lifefilters/medit.hpp>

/*
  
\author M.A. Fernandez
\date 18/10/2005

\brief  This programs couples the Navier-Stokes and (linear) Elastodynamic equations using the 
semi-implicit scheme introduced in  INRIA Research Report 5700, "A projection semi-implicit scheme for the coupling of 
an elastic structure with an incompressible fluid", by M.A. Fernández, J.-F. Gebeau and C. Grandmont.
Tested with SD stabilization (advective-diffusive sub-step) P1/P1 and Q1/Q1.

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
  
  
  //
  // Boundary conditions
  //
  BCHandler bcHu;  
  BCHandler bcHp; 
  BCHandler bcHd; 
  BCHandler bcHmesh;
  BCHandler bcHdp;
  BCHandler bcHdz;
  
  // pressure stress: computed not as a residual (implicitly coupledto the structre)
  //
  BCHandler sum_p;
  BCHandler sum_dp; // its derivative (for Newton on implicit sub-step)
  
  //========================================================================================
  // FLUID AND SOLID SOLVERS
  //========================================================================================
  //
  // The NavierStokes ALE solver
  //
  
  NavierStokesAleSolverCT< RegionMesh3D_ALE<LinearTetra> > fluid(data_file,
								 feTetraP1,
								 quadRuleTetra4pt,
								 quadRuleTria3pt, 
								 bcHu,
								 bcHp,
								 bcHmesh);
  
  // The structural solver
  //
  VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> > solid(data_file, 
							      feTetraP1, 
							      quadRuleTetra4pt,
							      quadRuleTria3pt, 
							      bcHd);
  
  
  //========================================================================================
  //  COUPLED FSI OPERATOR
  //========================================================================================
  //
  //
  
  
  operFS<  NavierStokesAleSolverCT< RegionMesh3D_ALE<LinearTetra> > ,  
    VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >  > 
    oper(fluid, solid, bcHdp, sum_p, sum_dp, bcHdz);
  
  
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
  // Passing data from the fluid to the structure
  //
  boost::shared_ptr<DofInterface3Dto3D> dofFluidToStructure( new DofInterface3Dto3D(feTetraP1, solid.dDof(), 
										      feTetraP1, fluid.uDof()) );
  dofFluidToStructure->update(solid.mesh(), 1, fluid.mesh(), 1, 0.0);
  
  
  // Passing data from structure to the fluid
  //
  boost::shared_ptr<DofInterface3Dto3D> dofStructureToFluid( new  DofInterface3Dto3D(feTetraP1, fluid.uDof(),
										     feTetraP1, solid.dDof()) );
  dofStructureToFluid->update(fluid.mesh(), 1, solid.mesh(), 1, 0.0);
  
  
  //========================================================================================
  //  BOUNDARY CONDITIONS
  //========================================================================================
  //
  
  
  BCFunctionBase bcf(fZero);
  BCFunctionBase p_in(g);
   
  BCVectorInterface displ(solid.disp(), dim_solid, dofStructureToFluid ); 
  BCVectorInterface g_wall(fluid.residual(), dim_fluid, dofFluidToStructure );
  
  BCVector u_wall(fluid.w(), dim_fluid);
  BCVectorInterface a_wall(oper.a(), dim_solid, dofStructureToFluid, 2); 
  BCVector p_wall(fluid.p(), dim_fluid, 1); 
 
  BCVectorInterface da_wall(oper.da(), dim_solid, dofStructureToFluid, 2); 
  BCVector dp_wall( fluid.dp(), dim_fluid, 1);
  
  sum_p.addBC("Interface" , 1, Natural, Full,  p_wall, 3); 
  sum_dp.addBC("Interface", 1, Natural, Full, dp_wall, 3); 
  
  
  // Boundary conditions for the harmonic extension of the
  // interface solid displacement
  bcHmesh.addBC("Interface", 1, Essential, Full, displ, 3);
  bcHmesh.addBC("Top",       3, Essential, Full, bcf,   3);
  bcHmesh.addBC("Base",      2, Essential, Full, bcf,   3);
  bcHmesh.addBC("Edges",    20, Essential, Full, bcf,   3);
  
  // Boundary conditions for the fluid velocity
  bcHu.addBC("Wall",        1, Essential, Full, u_wall, 3);    
  bcHu.addBC("Wall_Edges", 20, Essential, Full, u_wall, 3);
  bcHu.addBC("Inlet",       2, Natural,   Full, bcf   , 3);
  bcHu.addBC("Outlet",      3, Natural,   Full, bcf   , 3);

  // Boundary conditions for the pressure
  bcHp.addBC("Wall",        1, Natural  , Scalar, a_wall);    
  bcHp.addBC("Wall_Edges", 20, Essential, Scalar, p_in);
  bcHp.addBC("Inlet",       2, Essential, Scalar, p_in);
  bcHp.addBC("Outlet",      3, Essential, Scalar, bcf);
  
  
  // Boundary conditions for the solid displacement
  bcHd.addBC("Interface", 1, Natural, Full, g_wall, 3);
  bcHd.addBC("Top",       3, Essential, Full, bcf,  3);
  bcHd.addBC("Base",      2, Essential, Full, bcf,  3);

  // Boundary conditions for dp
  bcHdp.addBC("Wall",        1, Natural  , Scalar, da_wall);    
  bcHdp.addBC("Wall_Edges", 20, Essential, Scalar, bcf);
  bcHdp.addBC("Inlet",       2, Essential, Scalar, bcf);
  bcHdp.addBC("Outlet",      3, Essential, Scalar, bcf);

  // Boundary conditions for dz
  bcHdz.addBC("Interface", 1, Natural,   Full, g_wall, 3);
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
  
  Vector disp(3*dim_solid),  disp_1(3*dim_solid), 
    velo(3*dim_solid), velo_1(3*dim_solid); 
  
  disp = ZeroVector( disp.size() );
  disp_1 = ZeroVector( disp_1.size() );
  velo = ZeroVector( velo.size() );
  velo_1 = ZeroVector( velo_1.size() );

  ofstream out_iter("iter");
  ofstream out_res("res");
  
  // Temporal loop
  //
  for (Real time=dt; time <= T; time+=dt) 
    {
      
      fluid.timeAdvance( f, time);
      solid.timeAdvance(f,time);
      
      oper.setTime(time);

      //
      // EXPLICIT SUB-STEP: ALE reaction-advection-diffusion
      //
      
      fluid.updateMesh(time);
      fluid.iterateVelocity(time);
      
      
      //
      // IMPLCIT SUB-STEP: pressure + solid
      //
      Real dti=1.0/dt;
      
      velo = ( disp - disp_1 )*dti;	
      disp_1 = disp;
      
      disp = disp + dt*(1.5*velo - 0.5*velo_1);
      velo_1 = velo;
      
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
