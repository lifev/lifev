#include "lifeV.hpp"
#include "NavierStokesSolverPC.hpp"
#include "chrono.hpp" 
#include "ud_functions.hpp"
#include "GetPot.hpp"
#include <string> 

#undef MESH_INRIA
#define MESH_MOX
#undef FULLY_DIRICHLET

int main(int argc, char** argv)
{
  // Reading from data file 
  // 
  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);
  char name_wss[8];
  char name_sol[5];

  // Boundary conditions definition
  // 
  BCFunction_Base uin(u_in);
  BCFunction_Base uin_d(u_w);
  BCFunction_Base uw(u_w);
  BCFunction_Base uout(u_out);
  BCFunction_Base uoutd(u_w);
  
  BC_Handler BCh_u(7); // We impose 7 boundary conditions

  // mode "Full" involves all components, which are 3 for the velocity.
  BCh_u.addBC("Inlet",   10, Essential, Full, uin, 3);
  BCh_u.addBC("Inletd",  11, Essential, Full, uin_d, 3);
  BCh_u.addBC("Wall",    50, Essential, Full, uw, 3);
  BCh_u.addBC("Outlet1", 20, Natural,   Full, uout, 3);
  BCh_u.addBC("Outlet1d",21, Essential, Full, uoutd, 3);
  BCh_u.addBC("Outlet2", 30, Natural,   Full, uout, 3);
  BCh_u.addBC("Outlet2d",31, Essential, Full, uoutd, 3);
    
  //  BCFunction_Base u_wall(u1);
  //  BCFunction_Base in_flow(u0);
  //  BCFunction_Base out_flow(u2);
  //  BC_Handler BCh_u(3);  
  //  BCh_u.addBC("Wall",   20, Essential, Full, u_wall,  3);
  //  BCh_u.addBC("InFlow", 10, Essential,   Full, in_flow, 3);
  //  BCh_u.addBC("OutFlow", 50, Natural,   Full, out_flow, 3);
 

  // Navier-Stokes Solver
  //
  NavierStokesSolverPC< RegionMesh3D<LinearTetra> > ns(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt, 
						       quadRuleTria3pt, quadRuleTetra64pt, quadRuleTria3pt, BCh_u);
  ns.showMe();

  //  PhysVectUnknown<Vector> r(ns.uDof().numTotalDof());

  //ns.post_proc().show_bdLtoG();

  ns.post_proc_set_area();
  ns.post_proc_set_normal();
  ns.post_proc_set_phi();
 
//  ns.post_proc().show_area();
//  ns.post_proc().show_normal();
//  ns.post_proc().show_phi();
  
  // Initialization
  //
  Real dt = ns.timestep();
  Real T  = ns.endtime();
  ns.initialize(u0,p0,0.0,dt);  

  UInt index=0;  
  int nn;
    nn = sprintf(name_sol,"poi%d",index);
    ns.dx_write_sol(name_sol,"P1bubble","P1");
  // Temporal loop
  //
  for (Real time=dt ; time <= T; time+=dt) {
    index++;  
    ns.timeAdvance(f,time);
    ns.iterate(time);  
    ns.postProcessPressure();
    //    PhysVectUnknown<Vector> r(ns.residual());
    //    cout << "Residual" << endl;
    //    for (UInt index=0;index<r.size();index++) cout << index << " " << r.vec()(index) << endl;
    nn = sprintf(name_wss,"wss%d.dx",index);
    nn = sprintf(name_sol,"poi%d",index);
    ns.dx_write_sol(name_sol,"P1bubble","P1");
    ns.ShearStressCompute(name_wss,"P1bubble");
  }  
   
  return 0;
}
