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
  Real lambda0,Q,Qn2,Qn1,Qn;
  Real lambda;
  Real r0,y;
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

  // SECOND NS. I resolve once 
  // Boundary conditions definition
  //
  BCFunction_Base in_flow2(u2);
  BC_Handler BCh_u2(5);
  BCh_u2.addBC("Wall",   2, Essential, Full, u_wall,  3);
  BCh_u2.addBC("Wall-inflow",   4, Essential, Full, u_wall,  3);
  BCh_u2.addBC("Wall-outflow",   5, Essential, Full, u_wall,  3);
  BCh_u2.addBC("InFlow2", 1, Natural,   Full, in_flow2, 3);
  BCh_u2.addBC("OutFlow", 3, Natural,   Full, out_flow, 3);
   
  // Navier-Stokes Solver
  //
  NavierStokesSolverPC< RegionMesh3D<LinearTetra> > ns2(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt, 
						   quadRuleTria3pt, quadRuleTetra5pt, quadRuleTria3pt, BCh_u2);
  ns2.showMe();

  // Initialization
  //
  Real dt = ns2.timestep();  
  Real startT = ns2.inittime();
  Real T  = ns2.endtime();

  if(startT > 0.0){
     cout << "initialize velocity and pressure with data from file" << std::endl;
     ostringstream indexin;
     string vinname, cinname;
     indexin << (startT*100);
     vinname = "fluid.res"+indexin.str();
     ns2.initialize(vinname);}
  else{
     cout << "initialize velocity and pressure with u0 and p0" << std::endl;	
     ns2.initialize(u02,p0,0.0,dt);
  }

  // I solve only one time step because is always the same
  //
  for (Real time=startT+dt ; time <= startT+dt; time+=dt) {
    cout << endl;
    cout << "start NS2" << endl;
    ns2.timeAdvance(f,time);
    ns2.iterate(time); 
  }
 
  Qn2=ns2.flux(1); // compute the flux of NS2    
  cout << endl;
  cout << "end NS2" << endl;

  /////////
  // GMRes 
  /////////

  // Definitions to construct the vector lambda and the quantities for the exchange of informations between NS1 and NS3
  //
  PhysVectUnknown<Vector> u3(ns2.uDof().numTotalDof());
  UInt dim_lambda = ns2.uDof().numTotalDof(); 
  Vector vec_lambda(dim_lambda);
  BCVector bcvec(vec_lambda,dim_lambda,1);
  
  // Boundary conditions definition FIRST NS
  //
  BC_Handler BCh_u1(5);
  BCh_u1.addBC("Wall",   2, Essential, Full, u_wall,  3);
  BCh_u1.addBC("Wall-inflow",   4, Essential, Full, u_wall,  3);
  BCh_u1.addBC("Wall-outflow",   5, Essential, Full, u_wall,  3);
  BCh_u1.addBC("OutFlow", 3, Natural,   Full, out_flow, 3);
  BCh_u1.addBC("InFlow1", 1, Natural,   Full, bcvec, 3); 

  // Definition of the class for the FIRST NS
  //
  NavierStokesSolverPC< RegionMesh3D<LinearTetra> > ns1(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt, 
				       quadRuleTria3pt, quadRuleTetra5pt, quadRuleTria3pt, BCh_u1);
  ns1.showMe(); 

  // Initialization for FIRST NS
  //
  if(startT > 0.0){
     cout << "initialize velocity and pressure with data from file" << std::endl;
     ostringstream indexin;
     string vinname, cinname;
     indexin << (startT*100);
     vinname = "fluid.res"+indexin.str();
     ns1.initialize(vinname);}
  else{
     cout << "initialize velocity and pressure with u0 and p0" << std::endl;	
     ns1.initialize(u0,p0,0.0,dt);
  }
  
  Real init_flux1=ns1.flux(1);
  cout << endl;
  cout << "flusso iniziale NS1" << " " << init_flux1 << endl; 
  
  // Boundary conditions definition THIRD NS
  //
  BC_Handler BCh_u3(5);
  BCh_u3.addBC("Wall",   2, Essential, Full, u_wall,  3);
  BCh_u3.addBC("Wall-inflow",   4, Essential, Full, u_wall,  3);
  BCh_u3.addBC("Wall-outflow",   5, Essential, Full, u_wall,  3);
  BCh_u3.addBC("OutFlow", 3, Natural,   Full, out_flow, 3);
  BCh_u3.addBC("InFlow3", 1, Natural,   Full, bcvec, 3);

  // Definition of the class for the THIRD NS
  //
  NavierStokesSolverPC< RegionMesh3D<LinearTetra> > ns3(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt, 					       quadRuleTria3pt, quadRuleTetra5pt, quadRuleTria3pt, BCh_u3);
  ns3.showMe();

  // Initialization for THIRD NS
  //
  if(startT > 0.0){
     cout << "initialize velocity and pressure with data from file" << std::endl;
     ostringstream indexin;
     string vinname, cinname;
     indexin << (startT*100);
     vinname = "fluid.res"+indexin.str();
     ns3.initialize(vinname);}
  else{
     cout << "initialize velocity and pressure with u0 and p0" << std::endl;	
     ns3.initialize(u0,p0,0.0,dt);
  }

  Real init_flux3=ns3.flux(1);
  cout << endl;
  cout << "flusso iniziale NS3" << " " << init_flux3 << endl;
  
  //ofstream outfile("flusso.txt"); 

  // Temporal loop
  //
  for (Real time=startT+dt ; time <= T; time+=dt) {
 
    // Assigment of the flux and initialize lambda
    //
    Q=-1; 
    //Q=-1*cos(2*pi*time); // pulsatile flow
    //Q=-0.2*ns2.compute_flux(time); //fisiological flux
    lambda0=-Q;

    //outfile << Q << endl;

    // FIRST NS 
    // BC changing in time for first NS
    //
    for (UInt i=0 ; i < dim_lambda; ++i) {
      vec_lambda[i]=-lambda0;
    }
    
    // Navier-Stokes Solver
    //
    cout << endl;
    cout << "start NS1 time" << " " << time << endl;
    if (time==startT+dt) {
      ns1.timeAdvance(f,time);
    }
    else {
      ns1.timeAdvance(f,time,u3); // I pass the velocity u3 as the one at the previous time step 
    }
    ns1.iterate(time); 
    //
    // end FIRST NS
    
    Qn1=ns1.flux(1); //compute the flux of NS1

    // compute the variables to update lambda
    r0=Qn1-Q;
    if(r0>0){
      absr0=r0;
      v=1;
    }
    else{      
      absr0=-r0;
      v=-1;
      }
    y=absr0/Qn2;
    lambda=lambda0+v*y;    

    // THIRD NS
    //BC changing in time for THIRD NS
    //
    for (UInt i=0 ; i < dim_lambda; ++i) {
      vec_lambda[i]=-lambda;
    }
    
    // Navier-Stokes Solver
    //
    cout << "start NS3 time" << " " << time << endl;
    ns3.timeAdvance(f,time);
    u3=ns3.iterate(time,1); //return the velocity 

    Qn=ns3.flux(1); //compute the flux of NS3: the definitive one
    cout << "imposed flux" << " " << Q << endl;
    cout << "numerical flux" << " " << Qn << endl;

// ************* saving result on file *****************************************
    ostringstream indexout;
    indexout << (time*100);
    string voutname;
    voutname = "fluid.res"+indexout.str();
    fstream Resfile(voutname.c_str(),ios::out | ios::binary);
    Resfile.write((char*)&ns3.u()(1),ns3.u().size()*sizeof(double));
    Resfile.write((char*)&ns3.p()(1),ns3.p().size()*sizeof(double));
    Resfile.close();

    ns3.postProcess();
  }
  return 0;
}
