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
