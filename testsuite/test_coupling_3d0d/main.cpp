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
#include "zeroDModelSolver.hpp"



int main(int argc, char** argv)
{
  using namespace LifeV;

  // Reading from data file
  //
  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);

  BCFunctionBase u_wall(u1);
  BCFunctionBase out_flow(u1);
  BCHandler BCh_u(5);


  // Navier-Stokes Solver
  //
  typedef NavierStokesSolverPC< RegionMesh3D<LinearTetra> > NS;
  NS ns(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt,
	quadRuleTria3pt, quadRuleTetra5pt, quadRuleTria3pt, BCh_u);
  ns.showMe();

  UInt dim_fluid = ns.uDof().numTotalDof();
  Vector vec_press(dim_fluid);
  BCVector bcvec(vec_press,dim_fluid,1);

  BCh_u.addBC("Wall",   2, Essential, Full, u_wall,  3);
  BCh_u.addBC("Wall-inflow",   4, Essential, Full, u_wall,  3);
  BCh_u.addBC("Wall-outflow",   5, Essential, Full, u_wall,  3);
  BCh_u.addBC("InFlow", 1, Natural,  Full, bcvec , 3);
  BCh_u.addBC("OutFlow", 3, Natural,  Full, out_flow, 3);

  // Initialization
  //

  std::ofstream outfile;

  Real dt = ns.timestep();
  Real startT = ns.inittime();
  Real T  = ns.endtime();

  outfile.open("res_Q.m", std::ios::app);
  outfile << "xx=[0:" << dt << ":" << T << "]; " << std::endl;
  outfile.close();

  outfile.open("res_DP.m", std::ios::app);
  outfile << "xx=[0:" << dt << ":" << T << "]; " << std::endl;
  outfile.close();


  if(startT > 0.0){
     std::cout << "initialize velocity and pressure with data from file" << std::endl;
     std::ostringstream indexin;
     std::string vinname, cinname;
     indexin << (startT*100);
     vinname = "fluid.res"+indexin.str();
     ns.initialize(vinname);}
  else{
     std::cout << "initialize velocity and pressure with u0 and p0" << std::endl;
     ns.initialize(u0,p0,0.0,dt);
  }

  getPressureFromQ< NS > press(ns, data_file, 5, 0.5, 0);

  // Temporal loop
  //
  for (Real time=startT+dt ; time <= T; time+=dt) {

    vec_press = -press.pression(time);

    ns.timeAdvance(f,time);
    ns.iterate(time);

// ************* saving result on file *****************************************
    std::ostringstream indexout;
    indexout << (time*1000);
    std::string voutname;
    voutname = "fluid.res"+indexout.str();
    std::fstream Resfile(voutname.c_str(),std::ios::out | std::ios::binary);
    Resfile.write((char*)&ns.u()(1),ns.u().size()*sizeof(double));
    Resfile.write((char*)&ns.p()(1),ns.p().size()*sizeof(double));
    Resfile.close();



    ns.postProcess();
  }

  outfile.open("res_Q.m", std::ios::app);
  outfile << "    ]; " << std::endl;
  outfile << "figure;" << std::endl;
  outfile << "plot(xx,Q);" << std::endl;
  outfile.close();

  outfile.open("res_DP.m", std::ios::app);
  outfile << "     ]; " << std::endl;
  outfile << "figure;" << std::endl;
  outfile << "plot(xx,DP);" << std::endl;
  outfile.close();


  return 0;
}
