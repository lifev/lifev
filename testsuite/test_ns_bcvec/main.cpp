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
#include "NavierStokesSolverPC.hpp"
#include "chrono.hpp"
#include "ud_functions.hpp"
#include "GetPot.hpp"



int main(int argc, char** argv)
{
  using namespace LifeV;
  using namespace std;
  // Reading from data file
  //
  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);

  BCHandler BCh_u(5);

  // Navier-Stokes Solver
  //
  //Classical grad-grad version: put nothing or 0 in last argument
  //
  //NavierStokesSolverPC< RegionMesh3D<LinearTetra> > ns(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt,
  //				       quadRuleTria3pt, quadRuleTetra5pt, quadRuleTria3pt, BCh_u );
  //
  //
  // Curl-curl version: put 1 in last argument
  //
  NavierStokesSolverPC< RegionMesh3D<LinearTetra> > ns(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt,
						       quadRuleTria3pt, quadRuleTetra5pt, quadRuleTria3pt, BCh_u, 1 );
  ns.showMe();

  // boundary conditions definition
  //
  BCFunctionBase u_wall(u1);
  BCFunctionBase out_flow(u1);

  // Definition of the vector to impose as BC at inflow
  Real lambda=-2;
  UInt dim_lambda = ns.uDof().numTotalDof();
  Vector vec_lambda(dim_lambda);
  BCVector bcvec(vec_lambda,dim_lambda,1);

  for (UInt i=0 ; i < dim_lambda; ++i) {
    vec_lambda[i]=-lambda;
  }

  BCh_u.addBC("Wall",   2, Essential, Full, u_wall,  3);
  BCh_u.addBC("Wall-inflow",   4, Essential, Full, u_wall,  3);
  BCh_u.addBC("Wall-outflow",   5, Essential, Full, u_wall,  3);
  BCh_u.addBC("InFlow", 1, Natural,  Full, bcvec , 3);
  BCh_u.addBC("OutFlow", 3, Natural,  Full, out_flow, 3);


  // Initialization
  //
  Real dt = ns.timestep();
  Real startT = ns.inittime();
  Real T  = ns.endtime();

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


  Real p;

  // Temporal loop
  //
  for (Real time=startT+dt ; time <= T; time+=dt) {

    ns.timeAdvance(f,time);
    ns.iterate(time);

    // Calculation of the flux at inflow!
    cout << "flux  " << ns.flux(1);


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
