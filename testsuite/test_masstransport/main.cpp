// Main programm for coupled calculation of fluid-dynamics
// and mass transport in the arterial lumen  

// author:M. Prosi                                march/04
#include "lifeV.hpp"
#include "NavierStokesSolverPC.hpp"
#include "convDiffReactSolverPC.hpp" 
#include "chrono.hpp"
#include "ud_functions.hpp"
#include "GetPot.hpp"
#include "ensight7-write.hpp"

int main(int argc, char** argv)
{

 // ********** Reading from data file ******************************************

  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);
 

 // ********* Boundary conditions definition for fluid *************************
  BCFunction_Base u_wall(u1);
  BCFunction_Base u_inflow(u2);
  BC_Handler BCh_u(4); 

  BCh_u.addBC("Wall",   3, Essential, Full, u_wall,  3);  // non-permeable Wall - no-slip condition
  BCh_u.addBC("Wall-inflow",   4, Essential, Full, u_wall,  3);
  BCh_u.addBC("Wall-outflow",   5, Essential, Full, u_wall,  3);

  BCh_u.addBC("InFlow", 1, Essential,   Full, u_inflow, 3); // Velocity profile 

// ********** Boundary conditions definitions for mass transport ***************
  BCFunction_Base c_inflow(c1);
//  BCFunction_Mixte c_wall(alpha,beta);   // Permeability boundary condition
  BCFunction_Base c_wall(cw);              // Concentration boundary condition
  BC_Handler BCh_c(4);

//  BCh_c.addBC("C-Wall",   3, Mixte, Scalar, c_wall);  // Permeability boundary condition
//  BCh_c.addBC("C-Wall-inflow",   4, Mixte, Scalar, c_wall);
//  BCh_c.addBC("C-Wall-outflow",   5, Mixte, Scalar, c_wall);

  BCh_c.addBC("C-Wall",   3, Essential, Scalar, c_wall);  // Concentration boundary condition
  BCh_c.addBC("C-Wall-inflow",   4, Essential, Scalar, c_wall);
  BCh_c.addBC("C-Wall-outflow",   5, Essential, Scalar, c_wall);

  BCh_c.addBC("C-InFlow", 1, Essential,  Scalar, c_inflow); // Concentration profile 

// *********** Fluid class: ns *************************************************
  NavierStokesSolverPC< RegionMesh3D<LinearTetra> > 
	ns(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt, 
	quadRuleTria3pt, quadRuleTetra5pt, quadRuleTria3pt, BCh_u);
  ns.showMe();

// *********** Concentration class: cdr ****************************************
  ConvDiffReactSolverPC< RegionMesh3D<LinearTetra> > 
     cdr(data_file, feTetraP1, quadRuleTetra5pt, quadRuleTria3pt, BCh_c);
  cdr.showMe();

// *********** Initialization **************************************************
  Real dt = ns.timestep();
  Real T  = ns.endtime();
  Real startT = ns.inittime();

  if(startT > 0.0){
     std::cout << "initialize velocity and pressure with data from file" << std::endl;
     ostringstream indexin;
     std::string vinname, cinname;
     indexin << (startT*100);
     vinname = "fluid.res"+indexin.str();
     cinname = "concentration.res"+indexin.str();

     ns.initialize(vinname);
     cdr.initialize(cinname);
  }
  else{
     std::cout << "initialize velocity and pressure with u0 and p0" << std::endl;
     ns.initialize(u0,p0,0.0,dt);
     cdr.initialize(c0,0.0,dt);
  }

 // *** calculate the interpolation coordinates of the concentration grid points in the
 // *** velocity grid - for handling different grids for velocity and concenration
  ns.mesh().updateElementFaces(true);
  cdr.getcoord(ns.mesh(), ns.u(), BCh_u);

// ************ Temporal loop **************************************************
  for (Real time=startT+dt ; time <= T; time+=dt) {

     ns.timeAdvance(f,time);
     ns.iterate(time);

     cdr.timeAdvance(fc,time);

// ***** interpolate the velocity field on the concentration nodes ************
     cdr.getvel(ns.mesh(),ns.u(),BCh_u,time);

     cdr.iterate(time);

// ************* saving result on file *****************************************
    ostringstream indexout;
    indexout << (time*100);
    std::string voutname;
    voutname = "fluid.res"+indexout.str();
    std::fstream Resfile(voutname.c_str(),ios::out | ios::binary);
    Resfile.write((char*)&ns.u()(1),ns.u().size()*sizeof(double));
    Resfile.write((char*)&ns.p()(1),ns.p().size()*sizeof(double));
    Resfile.close();
    voutname = "concentration.res"+indexout.str();
    std::fstream Resfilec(voutname.c_str(),ios::out | ios::binary);
    Resfilec.write((char*)&cdr.c()(1),cdr.c().size()*sizeof(double));
    Resfilec.close();

// ************* creating Ensight output file **********************************
    outensight7Mesh3D(ns.mesh(), ns.u(), cdr.c(),time);

  }
  
  return 0;
}
