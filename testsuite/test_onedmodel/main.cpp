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
//! \author:Vincent Martin                                july/04
#include "lifeV.hpp"
#include "chrono.hpp"
#include "dataOneDModel.hpp"
#include "oneDModelSolver.hpp"
#include "ud_functions.hpp"
#include "GetPot.hpp"

int main(int argc, char** argv)
{

 //! ********** Reading from data file ******************************************

  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);

  OneDModelSolver onedm(data_file_name);
  onedm.showMeData();
  onedm.showMeHandler(cout, 6);

  return 0;
}
