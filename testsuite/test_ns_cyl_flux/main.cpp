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
#include "NavierStokesSolverPC.hpp"
#include "chrono.hpp"
#include "ud_functions.hpp"
#include "GetPot.hpp"

//////////////////////////////////////////////////////
// GMRES ALGHORITM FOR IMPOSED FLUX PROBLEM - 08/04 VC
//////////////////////////////////////////////////////

int main(int argc, char** argv)
{ 
  using namespace LifeV;
  using namespace std;

  // Definition of GMRes alghoritm variables and constants
  //
  Real lambda0,Q,Qno,Qn;
  Real lambda;
  Real r0,y,z;
  Real absr0;
  Real pi=3.14159265358979;
  int v; 
 
  // Reading from data file 
  //
  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);
 
  // BC common to all the NS
  //
  BCFunction_Base u_wall(u1);
  BCFunction_Base out_flow(u1);

  // Out of temporal loop NS. I resolve once 
  // Boundary conditions definition
  //
  BCFunction_Base in_flowo(uo);
  BC_Handler BCh_uo(5);
  BCh_uo.addBC("Wall",   2, Essential, Full, u_wall,  3);
  BCh_uo.addBC("Wall-inflow",   4, Essential, Full, u_wall,  3);
  BCh_uo.addBC("Wall-outflow",   5, Essential, Full, u_wall,  3);
  BCh_uo.addBC("InFlow2", 1, Natural,   Full, in_flowo, 3);
  BCh_uo.addBC("OutFlow", 3, Natural,   Full, out_flow, 3);
   
  // Navier-Stokes Solver
  //
  NavierStokesSolverPC< RegionMesh3D<LinearTetra> > nso(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra15pt, 
						   quadRuleTria3pt, quadRuleTetra5pt, quadRuleTria3pt, BCh_uo);
  nso.showMe();

  // Initialization
  //
  Real dt = nso.timestep();  
  Real startT = nso.inittime();
  Real T  = nso.endtime();

  if(startT > 0.0){
     cout << "initialize velocity and pressure with data from file" << std::endl;
     ostringstream indexin;
     string vinname, cinname;
     indexin << (startT*100);
     vinname = "fluid.res"+indexin.str();
     nso.initialize(vinname);}
  else{
     cout << "initialize velocity and pressure with u0 and p0" << std::endl;	
     nso.initialize(u0o,p0,0.0,dt);
  }

  // I solve only one time step because is always the same
  //
  for (Real time=startT+dt ; time <= startT+dt; time+=dt) {
    cout << endl;
    cout << "start NSo" << endl;
    nso.timeAdvance(f,time);
    nso.iterate(time); 
  }
 
  Qno=nso.flux(1); // compute the flux of NSo 
  cout << endl;
  cout << "end NSo" << endl;

  /////////
  // GMRes 
  /////////

  // Definitions to construct the vector lambda 
  //
  UInt dim_lambda = nso.uDof().numTotalDof(); 
  Vector vec_lambda(dim_lambda);
  BCVector bcvec(vec_lambda,dim_lambda,1);
  
  // Boundary conditions definition inner NS
  //
  BC_Handler BCh_u(5);
  BCh_u.addBC("Wall",   2, Essential, Full, u_wall,  3);
  BCh_u.addBC("Wall-inflow",   4, Essential, Full, u_wall,  3);
  BCh_u.addBC("Wall-outflow",   5, Essential, Full, u_wall,  3);
  BCh_u.addBC("OutFlow", 3, Natural,   Full, out_flow, 3);
  BCh_u.addBC("InFlow1", 1, Natural,   Full, bcvec, 3); 

  // Definition of the class
  //
  NavierStokesSolverPC< RegionMesh3D<LinearTetra> > ns(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra15pt, 
				       quadRuleTria3pt, quadRuleTetra5pt, quadRuleTria3pt, BCh_u);
  ns.showMe(); 

  // Initialization
  //
  if(startT > 0.0){
     cout << "initialize velocity and pressure with data from file" << std::endl;
     ostringstream indexin;
     string vinname, cinname;
     indexin << (startT*100);
     vinname = "fluid.res"+indexin.str();
     ns.initialize(vinname);}
  else{
     cout << "initialize velocity and pressure with u0 and p0" << std::endl;	
     ns.initialize(u0,p0,0.0,dt);
  }
  
  Qn=ns.flux(1);
  cout << endl;
  cout << "initial flux NS" << " " << Qn << endl; 
  
  //ofstream outfile("flusso.txt"); 

  // Temporal loop
  //
  for (Real time=startT+dt ; time <= T; time+=dt) {
 
    // Assigment of the flux and initialize lambda
    //
    Q=-1; 
    //Q=-1*cos(2*pi*time); // pulsatile flow
    //Q=-0.2*ns.compute_flux(time); //physiological flux
    lambda0=-Q;

    //outfile << Q << endl;

    // FIRST NS 
    // BC changing in time for NS
    //
    for (UInt i=0 ; i < dim_lambda; ++i) {
      vec_lambda[i]=-lambda0;
    }
    
    // Navier-Stokes Solver
    //
    cout << endl;
    cout << "start NS time" << " " << time << endl;
    ns.timeAdvance(f,time);
    ns.iterate(time); 
    //
    // end NS
    
    Qn=ns.flux(1); //compute the flux of NS

    // compute the variables to update lambda
    r0=Qn-Q;
    if(r0>0){
      absr0=r0;
      v=1;
    }
    else{      
      absr0=-r0;
      v=-1;
      }
    y=absr0/Qno;
    z=v*y;
    lambda=lambda0+z;    
 
    // update the velocity and the pressure
    ns.u()=ns.u()-z*nso.u();    
    ns.p()=ns.p()-z*nso.p();

    Qn=ns.flux(1); //compute the flux of NS: the definitive one
    cout << "imposed flux" << " " << Q << endl;
    cout << "numerical flux" << " " << Qn << endl;

// ************* saving result on file *****************************************
    ostringstream indexout;
    indexout << (time*100);
    string voutname;
    voutname = "fluid.res"+indexout.str();
    fstream Resfile(voutname.c_str(),ios::out | ios::binary);
    Resfile.write((char*)&ns.u()(1),ns.u().size()*sizeof(double));
    Resfile.write((char*)&ns.p()(1),ns.p().size()*sizeof(double));
    Resfile.close();

    ns.postProcess();
  }
  return 0;
}
