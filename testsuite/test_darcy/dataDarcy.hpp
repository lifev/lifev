#ifndef _DATADARCY_H_
#define _DATADARCY_H_
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "tab.hpp"

using namespace std;
class DataDarcy
{
public:
  //
  // Physics
  //
  int test_case;
  string mesh_file;
  int diffusion_type;
  double diffusion_coef;
  KNM<double> diffusion_tensor;
  //
  // Miscellaneous
  //
  string mesh_dir;
  string post_dir;
  int verbose;
  string post_proc_format;
  //
  DataDarcy(const GetPot& dfile);
  void dataDarcyShowMe(ostream& c=cout);
  void dataDarcyHelp(ostream& c=cout);
};
#endif
