#include <iostream>
#include "dataDarcy.hpp"

using namespace std;

//! constructor using a data file.
DataDarcy::DataDarcy(const GetPot& dfile):
  diffusion_tensor(3,3)
{
  // physics
  test_case = dfile("physics/test_case",1);
  mesh_file = dfile("physics/mesh_file","mesh.mesh");
  diffusion_type = dfile("physics/diffusion_type",0);
  diffusion_function = dfile("physics/diffusion_function",0);
  diffusion_scalar = dfile("physics/diffusion_scalar",1.);
  diffusion_tensor(0,0) = dfile("physics/diffusion_tensor",1.,0);
  diffusion_tensor(0,1) = dfile("physics/diffusion_tensor",0.,1);
  diffusion_tensor(0,2) = dfile("physics/diffusion_tensor",0.,2);
  diffusion_tensor(1,0) = dfile("physics/diffusion_tensor",0.,3);
  diffusion_tensor(1,1) = dfile("physics/diffusion_tensor",1.,4);
  diffusion_tensor(1,2) = dfile("physics/diffusion_tensor",0.,5);
  diffusion_tensor(2,0) = dfile("physics/diffusion_tensor",0.,6);
  diffusion_tensor(2,1) = dfile("physics/diffusion_tensor",0.,7);
  diffusion_tensor(2,2) = dfile("physics/diffusion_tensor",1.,8);
  // Linear Solver
  theLinearSolver = (LinearSolver) dfile("linearsolver/theLinearSolver",
                                         LinSlv_Aztec);
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
  c << "diffusion_scalar = " << diffusion_scalar << endl;
  c << "diffusion_tensor = " << diffusion_tensor << endl;
  // Linear Solver
  c << "\n*** Values for data [linear solver]\n\n";
  c << "theLinearSolver = " << theLinearSolver << endl;
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
  c << "test_case: a number indicating the test case" << endl;
  c << "mesh_file: the mesh file"<< endl;
  c << "diffusion_type: 0: scalar diffusion given by diffusion_scalar" << endl;
  c << "                1: tensor diffusion given diffusion_tensor" << endl;
  c << "                2: tensor diffusion for fibrous media" << endl;
  c << "diffusion_scalar: a scalar diffusion coefficient" << endl;
  c << "diffusion_tensor: the 9 entries of the diffusion tensor (by rows)" << endl;
  c << "diffusion_function: number of a user defined function (usrDiffusion.cc)" << endl;
  // Linear Solver
  c << "theLinearSolver: the linear solver you want to use:\n" 
    << "Aztec: " << LinSlv_Aztec <<  "\tUMFPack: " << LinSlv_UMFPack << endl;
  // miscellaneous
  c << "\n*** Help for data [miscellaneous]\n\n";
  c << "mesh_dir         : the directory where the mesh file is" << endl;;
  c << "post_dir         : the full postprocessing directory (including path)"
    << endl;
  c << "verbose          : to make the code verbose" << endl;
  c << "post_proc_format : postprocessing format (medit, vtk, ...)" << endl;
}

