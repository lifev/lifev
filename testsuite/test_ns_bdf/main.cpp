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
  BC_Handler BCh_u(2); 
  BCh_u.addBC("Wall",   2, Essential, Full, u_wall,  3);
  BCh_u.addBC("InFlow", 1, Natural,   Full, in_flow, 3);


  // Navier-Stokes Solver
  //
  NavierStokesSolverPC< RegionMesh3D<LinearTetra> > ns(data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt, 
						       quadRuleTria3pt, quadRuleTetra64pt, quadRuleTria3pt, BCh_u);
  ns.showMe();


  // Initialization
  //
  ns.initialize(u0);
  Real dt = ns.timestep();
  Real T  = ns.endtime();


  // Temporal loop
  //
  for (Real time=dt ; time <= T; time+=dt) {
    ns.timeAdvance(f,time);
    ns.iterate(time);  
    ns.postProcessPressure();
  }
  
  return 0;
}
