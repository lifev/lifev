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
#include <cstdlib>
#include "lifeV.hpp"
#include "norm.hpp"
#include "regionMesh3D_ALE.hpp"
#include "VenantKirchhofSolver.hpp"
#include "ud_functions.hpp"

/* 

   This program solves the linear elastodynamic equations (St Venant-Kirchhof material).
   The present test reproduces the free vibrations of a straight cylindrical vessel under
   a given initial perturbation. The strecture is fixed on its extremities.

*/


int main(int argc, char** argv)
{
  // Reading from data file 
  //
  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);
  
  // Number of boundary conditions for the velocity and mesh motion
  // 
  BC_Handler BCh(2); 

  // The linear Venant-Kirchhof solver
  //
  VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> > solid(data_file, feTetraP1, quadRuleTetra4pt, 
							  quadRuleTria3pt, BCh);
  solid.showMe();

  // Boundary conditions for the displacement
  //
  BCFunction_Base fixed(g1); 
  BCh.addBC("Base2 ", 2 , Essential, Full, fixed, 3); 
  BCh.addBC("Base3 ", 3 , Essential, Full, fixed, 3);

  // Temporal data and initial conditions
  //
  Real dt = solid.timestep();
  Real T  = solid.endtime();
  solid.initialize(d0,w0); // displacement and velocity
  
  // Temporal loop
  //
  for (Real time=dt; time <= T; time+=dt) {
    solid.timeAdvance(f,time); // Computes the rigth hand side
    solid.iterate(); // Computes the matrices and solves the system
    solid.postProcess(); // Post-presssing 
  }

  return 0;
}
