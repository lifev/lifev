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
#include <iostream>
#include "dataDarcy.hpp"

using namespace std;
DataDarcy::DataDarcy(const GetPot& dfile):
  diffusion_tensor(3,3)
{
  // physics
  test_case = dfile("physics/test_case",1);
  mesh_file = dfile("physics/mesh_file","mesh.me.hpp");
  diffusion_type = dfile("physics/diffusion_type",0);
  diffusion_coef = dfile("physics/diffusion_coef",1.);
  diffusion_tensor(0,0) = dfile("physics/diffusion_tensor",1.,0);
  diffusion_tensor(0,1) = dfile("physics/diffusion_tensor",0.,1);
  diffusion_tensor(0,2) = dfile("physics/diffusion_tensor",0.,2);
  diffusion_tensor(1,0) = dfile("physics/diffusion_tensor",0.,3);
  diffusion_tensor(1,1) = dfile("physics/diffusion_tensor",1.,4);
  diffusion_tensor(1,2) = dfile("physics/diffusion_tensor",0.,5);
  diffusion_tensor(2,0) = dfile("physics/diffusion_tensor",0.,6);
  diffusion_tensor(2,1) = dfile("physics/diffusion_tensor",0.,7);
  diffusion_tensor(2,2) = dfile("physics/diffusion_tensor",1.,8);
  // miscellaneous
  mesh_dir  = dfile("miscellaneous/mesh_dir","./");
  post_dir  = dfile("miscellaneous/post_dir","./");
  verbose   = dfile("miscellaneous/verbose",0);
  post_proc_format   = dfile("miscellaneous/post_proc_format","vtk");
}

void DataDarcy::dataDarcyShowMe(ostream& c)
{
  // physics
  c << "\n*** Values for data [physics]\n\n";
  c << "test_case = " << test_case << endl;
  c << "mesh_file = " << mesh_file << endl;
  c << "diffusion_type = " << diffusion_type << endl;
  c << "diffusion_coef = " << diffusion_coef << endl;
  c << "diffusion_tensor = " << diffusion_tensor << endl;
  // miscellaneous
  c << "\n*** Values for data [miscellaneous]\n\n";
  c << "mesh_dir         = " << mesh_dir << endl;
  c << "post_dir         = " << post_dir << endl;
  c << "verbose          = " << verbose << endl;
  c << "post_proc_format = " << post_proc_format << endl;
}

void DataDarcy::dataDarcyHelp(ostream& c)
{
  // physics
  c << "\n*** Help for data [physics]\n\n";
  c << "test_case: a number indicated the test case" << endl;
  c << "mesh_file: the mesh file"<< endl;
  c << "diffusion_type: 0: scalar diffusion, 1: tensor diffusion " << endl;
  c << "diffusion_coeff: a diffusion coefficient (the diffusion if diffusion_type = 0)" << endl;
  c << "diffusion_tensor: the 9 entries of the diffusion tensor (by rows)" << endl;
  // miscellaneous
  c << "\n*** Help for data [miscellaneous]\n\n";
  c << "mesh_dir         : the directory where the mesh file is" << endl;;
  c << "post_dir         : the directory where are saved the postprocessing files"
    << endl;
  c << "verbose          : to make the code verbose" << endl;
  c << "post_proc_format : postprocessing format (medit, vtk, ...)" << endl;
}

