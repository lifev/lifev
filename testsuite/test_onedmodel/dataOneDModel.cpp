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
/*!
  \file dataOneDModel.cpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief Implementation of the class for the data in oned model

*/
#include "dataOneDModel.hpp"

namespace LifeV
{
//! constructor using a data file.
DataOneDModel::DataOneDModel(const GetPot& dfile)
{
  //! Physics/BCs
  _M_test_case = dfile("physics/BCs/test_case",1);
  //! Physics/parameters
  _M_alphaCor = dfile("physics/parameters/alphaCor",1.);
  _M_beta0 = dfile("physics/parameters/beta0",1.);
  _M_beta1 = dfile("physics/parameters/beta1",0.5);
  _M_Kr = dfile("physics/parameters/Kr",1.);

  //! Time
  _M_time_beg = dfile("time/timebeg",0.0);
  _M_time_end = dfile("time/timeend",1.0);
  _M_time_step = dfile("time/timestep",0.1);
  //! Discretization
  _M_mesh_file = dfile("discretization/mesh_file","mesh.mesh");
  _M_mesh_dir  = dfile("discretization/mesh_dir","./");
  _M_x_left    = dfile("discretization/x_left",0.);
  _M_x_right   = dfile("discretization/x_right",1.);
  _M_nx        = dfile("discretization/nx",10);
  //! Miscellaneous
  _M_post_dir  = dfile("miscellaneous/post_dir","./");
  _M_verbose   = dfile("miscellaneous/verbose",0);
  _M_post_proc_format   = dfile("miscellaneous/post_proc_format","medit");
}

double DataOneDModel::timestep() const
{
  return _M_time_step;
}
double DataOneDModel::inittime() const
{
  return _M_time_beg;
}
double DataOneDModel::endtime() const
{
  return _M_time_end;
}

void DataOneDModel::showMeData(std::ostream& c) const
{
  c << "\n*** Values for data [physics]\n";
  // Physics/BCs
  c << "\t[physics/BCs]\n";
  c << "test_case = " << _M_test_case << "\n";
   //! Physics/parameters
  c << "\t[physics/parameters]\n";
  c << "alphaCor = " << _M_alphaCor << "\n";
  c << "beta0 = " << _M_beta0 << "\n";
  c << "beta1 = " << _M_beta1 << "\n";
  c << "Kr    = " << _M_Kr << "\n";

  //! Time
  c << "\n*** Values for data [time]\n";
  c << "time_beg = " << _M_time_beg << "\n";
  c << "time_end = " << _M_time_end << "\n";
  c << "time_step = " << _M_time_step << "\n";

  //! Discretization
  c << "\n*** Values for data [discretization]\n";
  c << "mesh_file = " << _M_mesh_file << "\n";
  c << "mesh_dir  = " << _M_mesh_dir << "\n";
  c << "x_left = " << _M_x_left << "\n";
  c << "x_right = " << _M_x_right << "\n";
  c << "nx = " << _M_nx << "\n";

  //! Miscellaneous
  c << "\n*** Values for data [miscellaneous]\n\n";
  c << "post_dir         = " << _M_post_dir << "\n";
  c << "verbose          = " << _M_verbose << "\n";
  c << "post_proc_format = " << _M_post_proc_format << std::endl;
}

}
