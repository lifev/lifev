#include <iostream>
#include "dataFhN.hpp"

namespace LifeV
{
  /*
    The following variables are static in order to make easier their
    passing to user functions (with a prototyp that cannot be changed)
  */
  int  DataFhN::test_case = 0;
  
  double DataFhN::stim_value_1 = 0.;
  double DataFhN::stim_radius_1 = 0.;
  double DataFhN::stim_start_1 = 0.;
  double DataFhN::stim_stop_1 = 0.;
  KN<double> DataFhN::stim_center_1(3);

  double DataFhN::stim_value_2 = 0.;
  double DataFhN::stim_radius_2 = 0.;
  double DataFhN::stim_start_2 = 0.;
  double DataFhN::stim_stop_2 = 0.;
  KN<double> DataFhN::stim_center_2(3);



  //! constructor using a data file.
  DataFhN::DataFhN(const GetPot& dfile)
  {
    // physics
    test_case = dfile("physics/test_case",1);
    mesh_file = dfile("physics/mesh_file","mesh.mesh");
    fhn_diff = dfile("physics/fhn_diff",1.);
    fhn_diff_fct = dfile("physics/fhn_diff_fct",0);
    fhn_f0 = dfile("physics/fhn_f0",1.);
    fhn_eps =  dfile("physics/fhn_eps",0.01);
    fhn_alpha   = dfile("physics/fhn_alpha",0.25);
    fhn_beta   = dfile("physics/fhn_beta",0.25);
    fhn_gamma   = dfile("physics/fhn_gamma",0.25);

    t_start_LV = dfile("physics/t_start_LV",0.);
    t_stop_LV = dfile("physics/t_stop_LV",0.);
    source_value_LV = dfile("physics/source_value_LV",0.);
    
    t_start_RV = dfile("physics/t_start_RV",0.);
    t_stop_RV = dfile("physics/t_stop_RV",0.);
    source_value_RV = dfile("physics/source_value_RV",0.);


    stim_start_1 = dfile("physics/stim_start_1",0.);
    stim_stop_1 = dfile("physics/stim_stop_1",0.);
    stim_value_1 = dfile("physics/stim_value_1",0.);
    DataFhN::stim_radius_1 = dfile("physics/stim_radius_1",0.);
    stim_center_1(0) =  dfile("physics/stim_center_1",0.,0);
    stim_center_1(1) =  dfile("physics/stim_center_1",0.,1);
    stim_center_1(2) =  dfile("physics/stim_center_1",0.,2);

    stim_start_2 = dfile("physics/stim_start_2",0.);
    stim_stop_2 = dfile("physics/stim_stop_2",0.);
    stim_value_2 = dfile("physics/stim_value_2",0.);
    stim_radius_2 = dfile("physics/stim_radius_2",0.);
    stim_center_2(0) =  dfile("physics/stim_center_2",0.,0);
    stim_center_2(1) =  dfile("physics/stim_center_2",0.,1);
    stim_center_2(2) =  dfile("physics/stim_center_2",0.,2);

    
    // miscellaneous
    mesh_dir  = dfile("miscellaneous/mesh_dir","./");
    post_dir  = dfile("miscellaneous/post_dir","./");
    verbose   = dfile("miscellaneous/verbose",0);
    post_proc_format   = dfile("miscellaneous/post_proc_format","vtk");
    nb_ref_time   = dfile("miscellaneous/nb_ref_time",10);
    ref_file   = dfile("miscellaneous/ref_file","ref.mat");
    store_file   = dfile("miscellaneous/store_file","sol.mat");
  }

void DataFhN::dataFhNShowMe(std::ostream& c)
{
  // physics
  c << "\n*** Values for data [physics]\n\n";
  c << "test_case = " << test_case << std::endl;
  c << "mesh_file = " << mesh_file << std::endl;
  c << "fhn_diff = " << fhn_diff << std::endl;
  c << "fhn_diff_fct = " << fhn_diff_fct << std::endl;
  c << "fhn_f0 = " << fhn_f0 << std::endl;
  c << "fhn_eps = " << fhn_eps << std::endl;
  c << "fhn_alpha = " << fhn_alpha << std::endl;
  c << "fhn_beta = " << fhn_beta << std::endl;
  c << "fhn_gamma = " << fhn_beta << std::endl;

  c << "t_start_LV = " << t_start_LV << std::endl;
  c << "t_stop_LV = " << t_stop_LV << std::endl;
  c << "source_value_LV = " << source_value_LV << std::endl;
  c << "stim_start_1 = " << stim_start_1 << std::endl;
  c << "stim_stop_1 = " << stim_stop_1 << std::endl;
  c << "stim_value_1 = " << stim_value_1 << std::endl;
  c << "stim_radius_1 = " << DataFhN::stim_radius_1 << std::endl;
  c << "stim_center_1 = "
    << stim_center_1(0) << ","
    << stim_center_1(1) << "," 
    << stim_center_1(2) << std::endl;
  
  c << "t_start_RV = " << t_start_RV << std::endl;
  c << "t_stop_RV = " << t_stop_RV << std::endl;
  c << "source_value_RV = " << source_value_RV << std::endl;
  c << "stim_start_2 = " << stim_start_2 << std::endl;
  c << "stim_stop_2 = " << stim_stop_2 << std::endl;
  c << "stim_value_2 = " << stim_value_2 << std::endl;
  c << "stim_radius_2 = " << stim_radius_2 << std::endl;
  c << "stim_center_2 = "
    << stim_center_2(0) << ","
    << stim_center_2(1) << "," 
    << stim_center_2(2) << std::endl;

  // miscellaneous
  c << "\n*** Values for data [miscellaneous]\n\n";
  c << "mesh_dir         = " << mesh_dir << std::endl;
  c << "post_dir         = " << post_dir << std::endl;
  c << "verbose          = " << verbose << std::endl;
  c << "post_proc_format = " << post_proc_format << std::endl;
  c << "nb_ref_time      = " << nb_ref_time << std::endl;
  c << "ref_file      = " << ref_file << std::endl;
  c << "store_file      = " << store_file << std::endl;
}

void DataFhN::dataFhNHelp(std::ostream& c)
{
  // physics
  c << "\n*** Help for data [physics]\n\n";
  c << "test_case: a number indicating the test case" << std::endl;
  c << "mesh_file: the mesh file"<< std::endl;
  c << "fhn_diff = the diffusion constant in the Fitzhugh-Nagumo equations"  << std::endl;
  c << "fhn_diff_fct = 0 : constant diffusion, 1 : space dependent diffusion"  << std::endl;
  c << "fhn_f0 = the constant before the cubic term in the Fitzhugh-Nagumo equations"
    << std::endl;
  c << "fhn_eps = the epsilon in the Fitzhugh-Nagumo equations" << std::endl;
  c << "fhn_alpha = the constant 'alpha' in the definition of the cubic"  << std::endl;
  c << "fhn_beta = the beta constant in the Fitzhugh-Nagumo equations"  << std::endl;
  c << "fhn_gamma = the gamma constant in the Fitzhugh-Nagumo equations"  << std::endl;
  c << "t_start_RV = starting time for the source term, in the RV" << std::endl;
  c << "t_stop_RV = stopping time for the source term, in the RV" << std::endl;
  c << "source_value_RV = value of the source term (inside a sphere)" << std::endl;


  c << "t_start_LV = starting time for the source term, in the LV" << std::endl;
  c << "t_stop_LV = stopping time for the source term, in the LV" << std::endl;
  c << "source_value_LV = value of the source term (inside a sphere)" << std::endl;

  c << "stim_start_1 = starting time for the stimulation 1" << std::endl;
  c << "stim_stop_1 = stopping time for the stimulation  2" << std::endl;
  c << "stim_value_1 = value of the stimulation 1 (inside a sphere)" << std::endl;
  c << "stim_radius_1 = radius of the sphere for the stimulation  1"
    << std::endl;
  c << "stim_center_1 = x,y,z of the center of the sphere for the stimulation 1" << std::endl;

  c << "stim_start_2 = starting time for the stimulation 2" << std::endl;
  c << "stim_stop_2 = stopping time for the stimulation  2" << std::endl;
  c << "stim_value_2 = value of the stimulation 2 (inside a sphere)" << std::endl;
  c << "stim_radius_2 = radius of the sphere for the stimulation  2"
    << std::endl;
  c << "stim_center_2 = x,y,z of the center of the sphere for the stimulation 2" << std::endl;

  c << "nb_ref_time = number of time step stored in/for the reference solution";
  c << "ref_file = file where a reference solution is load";
  c << "store_file = file where the current solution is stored";
  // miscellaneous
  c << "\n*** Help for data [miscellaneous]\n\n";
  c << "mesh_dir         : the directory where the mesh file is" << std::endl;;
  c << "post_dir         : the full postprocessing directory (including path)"
    << std::endl;
  c << "verbose          : to make the code verbose" << std::endl;
  c << "post_proc_format : postprocessing format (medit, vtk, ...)" << std::endl;
}
}
