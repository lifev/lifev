#define MESH_INRIA
#include "lifeV.hpp"
#include "NavierStokesAleSolverPC.hpp"
#include "VenantKirchhofSolver.hpp"
#include "ud_functions.hpp"
#include "picard.hpp"
#include "operFS.hpp"
#include "norm.hpp"
#include "dofInterface3Dto3D.hpp"


/*

   This programs couples the Navier-Stokes and (linear) Elastodynamic equations
   At each time step the resulting non-linear coupled problem is solved via
   Picard iterations with Aitken's acceleration. 

   The present test simulates the pressure wave propagation in a straight cylindrical vessel

*/


int main(int argc, char** argv)
{
  // Reading from data file 
  //
  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);

   
 
  // Number of boundary conditions for the fluid velocity, 
  // solid displacement, and fluid mesh motion
  // 
  BC_Handler BCh_u(3,0); 
  BC_Handler BCh_d(3,0);   
  BC_Handler BCh_mesh(4,1); 


  //========================================================================================
  // FLUID AND SOLID SOLVERS
  //========================================================================================
  //
  // The NavierStokes ALE solver
  //
  NavierStokesAleSolverPC< RegionMesh3D<LinearTetra> > fluid(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt, 
							     quadRuleTria3pt, quadRuleTetra64pt, quadRuleTria3pt, 
							     BCh_u,BCh_mesh);

  // The structural solver
  //
  VenantKirchhofSolver< RegionMesh3D<LinearTetra> > solid(data_file, feTetraP1, quadRuleTetra4pt, 
							  quadRuleTria3pt, BCh_d);

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
  BCVector_Interface g_wall(fluid.residual(), dim_fluid, dofFluidToStructure);

 
  // Passing data from structure to the fluid mesh: motion of the fluid domain
  //
  DofInterface3Dto3D dofStructureToFluidMesh(fluid.mesh().getRefFE(), fluid.dofMesh(), 
				       feTetraP1, solid.dDof()); 
  dofStructureToFluidMesh.update(fluid.mesh(), 1, solid.mesh(), 1, 0.0);                                                     
  BCVector_Interface displ(solid.d().vec(), dim_solid, dofStructureToFluidMesh);
 


  // Passing data from structure to the fluid: solid velocity at the interface velocity
  //
  DofInterface3Dto3D dofMeshToFluid(feTetraP1bubble, fluid.uDof(), feTetraP1bubble, fluid.uDof() );
  dofMeshToFluid.update(fluid.mesh(), 1, fluid.mesh(), 1, 0.0);                                                     
  BCVector_Interface u_wall(fluid.wInterpolated().vec(),fluid.uDof().numTotalDof(),dofMeshToFluid);


  //========================================================================================
  //  BOUNDARY CONDITIONS
  //========================================================================================
  //
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
  BCh_d.addBC("Interface", 1, Natural, Full, g_wall, 3);
  BCh_d.addBC("Top",       3, Essential, Full, bcf,  3);
  BCh_d.addBC("Base",      2, Essential, Full, bcf,  3);
 

  //========================================================================================
  //  COUPLED FSI OPERATOR
  //========================================================================================
  //
  operFS oper(fluid,solid);
 
  //========================================================================================
  //  TEMPORAL LOOP
  //========================================================================================
  //
  
  UInt maxpf = 100;
  Real dt = fluid.timestep();
  Real T  = fluid.endtime();
  fluid.initialize(u0);
  solid.initialize(d0,w0);
  Real omega=0, abstol=1.e-7, reltol=0.0;
  int status,maxiter;
 
  ofstream nout("num_iter");
  ASSERT(nout,"Error: Output file cannot be opened.");


  Vector disp(3*dim_solid);
  disp =0.0;
  Vector disp_old(3*dim_solid);
  disp_old =0.0;
  Vector dispStruct(3*dim_solid);
  dispStruct =0.0;
  Vector dispStruct_old(3*dim_solid);
  dispStruct_old =0.0;
  Vector velo(3*dim_solid);
  velo =0.0; 
  Vector velo_old(3*dim_solid);
  velo_old =0.0;
  Vector velo_1(3*dim_solid);
  velo_1 =0.0;

  for (Real time=dt; time <= T; time+=dt) {
    
    fluid.timeAdvance(f,time); 
    solid.timeAdvance(f,time);

    // Displacement prediction
    //
    disp = dispStruct + dt*(1.5*velo - 0.5*velo_1);
    
    velo_1 = velo;
    
    disp_old = disp;
    cout << "        norm(dispStruct  ) init = " << maxnorm(disp) << endl;
    cout << "        norm(velo  ) init = "       << maxnorm(velo) << endl;
    cout << "        norm(velo_1) init = "       << maxnorm(velo_1) << endl;
   
    maxiter = maxpf;

    // Picard-Aitken iterations
    //
    status = picard(&oper,maxnorm,dispStruct,dispStruct_old,velo,velo_old,
		    disp,disp_old,abstol,reltol,maxiter,1,omega);
    
    if(status == 1) {
      cout << "Inners iterations failed\n";
      exit(1);
    }  else {
      cout << "End of time "<< time << endl; 
      nout << time << "   " << maxiter << endl; 
      cout << "Number of inner iterations       : " << maxiter << endl;
      fluid.postProcess();  
      solid.postProcess(); 
    }
  }
  return 0;
}
   


