/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <iostream>
#include "dataTransient.hpp"

/*
  \author JFG
  \brief  poor data for transient dependent problem
  \TODO
  1) allow variable time steps 
  2) select a stopping test (based on either max_time_iter or max_time)
  3) a standard banner for new time step
  4) tolerance for steady state
*/

namespace LifeV
{
  
  using namespace std;
  
  DataTransient::DataTransient(const GetPot& dfile)
  {
    time_step = dfile("transient/time_step",0.1);
    max_time_iter = dfile("transient/max_time_iter",10);
    max_time      = dfile("transient/max_time",1.);
    post_proc_period      = dfile("transient/post_proc_period",1);
  }
  
  void DataTransient::dataTransientShowMe(ostream& c)
  {
    c << "time_step      = " << time_step << endl;
    c << "max_time_iter  = " << max_time_iter << endl;
    c << "max_time       = " << max_time << endl;
    c << "post_proc_period = " << post_proc_period << endl;
  }

  void DataTransient::dataTransientHelp(ostream& c)
  {
    c << "\n*** Help for data [transient]\n\n";
    c << "time_step        =  the time step" << endl;
    c << "max_time_iter    = maximum number of time iterations" << endl;
    c << "max_time         = maximum time" << endl;
    c << "post_proc_period  = periodicity of the post-processing (number of iterations)" << endl;
  }
  void DataTransient::timeBanner(int iter,double time)
  {
    cout << "======================================================================\n";
    cout << "Iter = " << iter << " Time = " << time << endl;
    cout << "======================================================================\n";
  }
  
}


