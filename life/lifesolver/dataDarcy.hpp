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
#ifndef _DATADARCY_H_
#define _DATADARCY_H_
#include <string>
#include <iostream>
#include "dataMesh.hpp"
#include "GetPot.hpp"
#include "tab.hpp"

namespace LifeV
{

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


template <typename Mesh>
class DataDarcy:
        public DataMesh<Mesh>
{
public:
    //
    // Physics
    //
    int test_case;
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
    std::string post_dir; //! full name (including path)
    int verbose;
    std::string post_proc_format;

    //! constructor using a data file.
    DataDarcy(const GetPot& dfile);
    void dataDarcyShowMe(std::ostream& c=std::cout);
    void dataDarcyHelp(std::ostream& c=std::cout);
};


// --------------------------------------------------
// IMPLEMENTATION
// --------------------------------------------------

//! constructor using a data file.
template <typename Mesh>
DataDarcy<Mesh>::DataDarcy(const GetPot& dfile):
    DataMesh<Mesh>( dfile, "darcy/discretization" ),
    diffusion_tensor(3,3)
{
    // physics
    test_case = dfile("darcy/physics/test_case",1);
    diffusion_type = dfile("darcy/physics/diffusion_type",0);
    diffusion_function = dfile("darcy/physics/diffusion_function",0);
    diffusion_scalar = dfile("darcy/physics/diffusion_scalar",1.);
    diffusion_tensor(0,0) = dfile("darcy/physics/diffusion_tensor",1.,0);
    diffusion_tensor(0,1) = dfile("darcy/physics/diffusion_tensor",0.,1);
    diffusion_tensor(0,2) = dfile("darcy/physics/diffusion_tensor",0.,2);
    diffusion_tensor(1,0) = dfile("darcy/physics/diffusion_tensor",0.,3);
    diffusion_tensor(1,1) = dfile("darcy/physics/diffusion_tensor",1.,4);
    diffusion_tensor(1,2) = dfile("darcy/physics/diffusion_tensor",0.,5);
    diffusion_tensor(2,0) = dfile("darcy/physics/diffusion_tensor",0.,6);
    diffusion_tensor(2,1) = dfile("darcy/physics/diffusion_tensor",0.,7);
    diffusion_tensor(2,2) = dfile("darcy/physics/diffusion_tensor",1.,8);
    // Linear Solver
    theLinearSolver = (LinearSolver) dfile("darcy/linearsolver/theLinearSolver",
                                           LinSlv_Aztec);
    // miscellaneous
    post_dir  = dfile("darcy/miscellaneous/post_dir","./");
    verbose   = dfile("darcy/miscellaneous/verbose",0);
    post_proc_format   = dfile("darcy/miscellaneous/post_proc_format","medit");
}

template <typename Mesh>
void DataDarcy<Mesh>::dataDarcyShowMe(std::ostream& c)
{
    // physics
    c << "\n*** Values for data [physics]\n\n";
    c << "test_case = " << test_case << "\n";
    c << "diffusion_type = " << diffusion_type << "\n";
    c << "diffusion_scalar = " << diffusion_scalar << "\n";
    c << "diffusion_tensor = " << diffusion_tensor << "\n";
    // Linear Solver
    c << "\n*** Values for data [linear solver]\n\n";
    c << "theLinearSolver = " << theLinearSolver << "\n";
    // miscellaneous
    c << "\n*** Values for data [miscellaneous]\n\n";
    c << "post_dir         = " << post_dir << "\n";
    c << "verbose          = " << verbose << "\n";
    c << "post_proc_format = " << post_proc_format << std::endl;
}

template <typename Mesh>
void DataDarcy<Mesh>::dataDarcyHelp(std::ostream& c)
{
    // physics
    c << "\n*** Help for data [physics]\n\n";
    c << "test_case: a number indicating the test case" << "\n";
    c << "diffusion_type: 0: scalar diffusion given by diffusion_scalar" << "\n";
    c << "                1: tensor diffusion given diffusion_tensor" << "\n";
    c << "                2: tensor diffusion for fibrous media" << "\n";
    c << "diffusion_scalar: a scalar diffusion coefficient" << "\n";
    c << "diffusion_tensor: the 9 entries of the diffusion tensor (by rows)" << "\n";
    c << "diffusion_function: number of a user defined function (usrDiffusion.cc)" << "\n";
    // Linear Solver
    c << "theLinearSolver: the linear solver you want to use:\n"
      << "Aztec: " << LinSlv_Aztec <<  "\tUMFPack: " << LinSlv_UMFPack << "\n";
    // miscellaneous
    c << "\n*** Help for data [miscellaneous]\n\n";
    c << "post_dir         : the full postprocessing directory (including path)"
      << "\n";
    c << "verbose          : to make the code verbose" << "\n";
    c << "post_proc_format : postprocessing format (medit, vtk, ...)" << std::endl;
}

}
#endif
