#define MESH_INRIA
#include "lifeV.hpp"
#include "NavierStokesSolverPC.hpp"
#include "chrono.hpp"
#include "ud_functions.hpp"
#include "GetPot.hpp"



int main(int argc, char** argv)
{
  // Reading from data file 
  //
  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);
 

  // Boundary conditions definition
  //
  BCFunction_Base u_wall(u1);
  BCFunction_Base in_flow(u2);
  BC_Handler BCh_u(4);
  BCh_u.addBC("Wall",   2, Essential, Full, u_wall,  3);
  BCh_u.addBC("Wall-inflow",   4, Essential, Full, u_wall,  3);
  BCh_u.addBC("Wall-outflow",   5, Essential, Full, u_wall,  3);
  BCh_u.addBC("InFlow", 1, Essential,   Full, in_flow, 3);
   


  // Navier-Stokes Solver
  //
  NavierStokesSolverPC< RegionMesh3D<LinearTetra> > ns(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt, 
						       quadRuleTria3pt, quadRuleTetra64pt, quadRuleTria3pt, BCh_u);
  ns.showMe();


  // Initialization
  //
  Real dt = ns.timestep();
  Real T  = ns.endtime();
  ns.initialize(u0,p0,0.0,dt);


  // Temporal loop
  //
  for (Real time=dt ; time <= T; time+=dt) {
    ns.timeAdvance(f,time);
    ns.iterate(time);  
    ns.postProcessPressure();
  }
  
  return 0;
}
