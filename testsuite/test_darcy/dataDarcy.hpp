// data darcy declaration -*- C++ -*-

#ifndef _DATADARCY_H_
#define _DATADARCY_H_
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "tab.hpp"

namespace LifeV
{
using namespace std;

//! add here any linear solver you like (Mumps...)
enum LinearSolver{LinSlv_Aztec, LinSlv_UMFPack};

/*! Type of BC on the interface (if any) of the subdomain
 in Input (In) or Output (Out).
 If p is the scalar variable and u = K grad p the velocity
 and n the unit normal
   -- Dirichlet  : p                    on Sigma
   -- Neumann    : (- u.n)              on Sigma
   -- Robin      : (- u.n) + alpha p    on Sigma
   -- MinusRobin : -(- u.n) + alpha p   on Sigma
*/
enum TypeOfBCInterfaceIn{DirichletIn, NeumannIn, RobinIn};
enum TypeOfBCInterfaceOut{DirichletOut, NeumannOut, RobinOut, MinusRobinOut};


class DataDarcy
{
public:
  //
  // Physics
  //
  int test_case;
  string mesh_file;
  int diffusion_type;
  int diffusion_function;
  double diffusion_scalar;
  KNM<double> diffusion_tensor;
  //
  // Linear Solver
  //
  LinearSolver theLinearSolver;
  //
  // Miscellaneous
  //
  string mesh_dir;
  string post_dir; //! full name (including path)
  int verbose;
  string post_proc_format;

  //! constructor using a data file.
  DataDarcy(const GetPot& dfile);
  void dataDarcyShowMe(ostream& c=cout);
  void dataDarcyHelp(ostream& c=cout);
};
}
#endif
