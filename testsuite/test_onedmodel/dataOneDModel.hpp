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
  \file dataOneDModel.hpp
  \author Vincent Martin
  \date 07/2004 
  \version 1.0

  \brief File containing a class for data in oned model 

*/
#ifndef _DATAONEDMODEL_H_
#define _DATAONEDMODEL_H_
#include <string>
#include <iostream>
#include "lifeV.hpp"
#include "GetPot.hpp"


/*! 
  \class DataOneDModel

  Base class which holds usual data for the Blood flow One Dimensional Model solvers

*/
class DataOneDModel 
{
public:

  //! Constructor
  DataOneDModel(const GetPot& dfile);
  
  //! Output 
  void showMeData(std::ostream& c=std::cout);

protected:
  //
  //! Physics
  //
  //! Physics/BCs (boundary conditions)
  int _M_test_case; //! test case
  //! Physics/parameters
  double _M_alphaCor; //! coriolis coeff
  double _M_beta0;  //! pressure law factor
  double _M_beta1;  //! pressure law exponent
  double _M_Kr;     //! friction parameter (source term)

  //
  //! Time
  //
  double _M_time_beg;  //! starting time
  double _M_time_end; //! finishing time
  double _M_time_step; //! time step
  //
  //! Discretization
  //
  std::string _M_mesh_file;
  std::string _M_mesh_dir;
  double _M_x_left;  //! left coordinate
  double _M_x_right; //! right coordinate
  int _M_nx; //! such that (nx+1)=number of elements (0,1,..,nx)
  //
  //! Miscellaneous
  //
  std::string _M_post_dir; //! full name (including path)
  int _M_verbose;
  std::string _M_post_proc_format;
};

#endif
