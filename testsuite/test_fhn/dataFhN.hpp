// data fhn declaration -*- C++ -*-

#ifndef _DATAFHN_H_
#define _DATAFHN_H_
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "tab.hpp"

namespace LifeV
{
  using namespace std;
  
  class DataFhN
  {
  public:
    //
    // Physics
    //
    /*    
	  du/dt - diff \lapl u = f0 u(1-u)(u-alpha) - v
	  dv/dt = eps (beta u - gamma v)
    */
    static int test_case;
    string mesh_file;
    double fhn_diff;
    double fhn_f0;
    double fhn_eps;
    double fhn_alpha;
    double fhn_beta;
    double fhn_gamma;
    double t_start_LV;
    double t_stop_LV;
    double source_value_LV;
    double t_start_RV;
    double t_stop_RV;
    double source_value_RV;
    static KN<double> stim_center_1;
    static double stim_radius_1;
    static double stim_start_1;
    static double stim_stop_1;
    static double stim_value_1;
    static KN<double> stim_center_2;
    static double stim_radius_2;
    static double stim_start_2;
    static double stim_stop_2;
    static double stim_value_2;
    int nb_ref_time;
    string ref_file;
    string store_file;
    int fhn_diff_fct; 
    //
    // Miscellaneous
    //
    string mesh_dir;
    string post_dir; //! full name (including path)
    int verbose;
    string post_proc_format;
    //! constructor using a data file.
    DataFhN(const GetPot& dfile);
    void dataFhNShowMe(ostream& c=cout);
    void dataFhNHelp(ostream& c=cout);
  };
}
#endif
