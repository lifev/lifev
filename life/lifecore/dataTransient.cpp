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
  
  DataTransient::DataTransient(const GetPot& dfile)
  {
    init_time = dfile("transient/init_time",0.);
    time_step = dfile("transient/time_step",0.1);
    max_time_iter = dfile("transient/max_time_iter",10);
    max_time      = dfile("transient/max_time",1.);
    post_proc_period      = dfile("transient/post_proc_period",1);
    adapt_period      = dfile("transient/adapt_period",1);
    init_data =  dfile("transient/init_data",0);
    init_file_name = dfile("transient/init_file_name","sol"); 
  }
  
  void DataTransient::dataTransientShowMe(std::ostream& c)
  {
    c << "time_step      = " << time_step << std::endl;
    c << "max_time_iter  = " << max_time_iter << std::endl;
    c << "init_time       = " << init_time << std::endl;
    c << "max_time       = " << max_time << std::endl;
    c << "post_proc_period = " << post_proc_period << std::endl;
    c << "adapt_period = " << adapt_period << std::endl;
    c << "init_data  = " << init_data << std::endl;
    c << "init_file_name  = " << init_file_name << std::endl;
  }


void DataTransient::dataTransientHelp( std::ostream& c )
{
    c << "\n*** Help for data [transient]\n\n";
    c << "time_step        =  the time step" << std::endl;
    c << "max_time_iter    = maximum number of time iterations" << std::endl;
    c << "init_time         = initial time" << std::endl;
    c << "max_time         = maximum time" << std::endl;
    c << "post_proc_period  = periodicity of the post-processing (number of iterations)" << std::endl;
    c << "adapt_period  = periodicity for the adaptation (number of iterations)" << std::endl;
    c << "init_data  = 0 if the initial solution is computed in the code, 1 if the initial solution is read on a file"  << std::endl;
    c << "init_file_name  = the name of the file where the initial solution is read"  << std::endl;
  }
  void DataTransient::timeBanner(int iter,double time)
  {
    std::cout << "======================================================================\n";
    std::cout << "Iter = " << iter << " Time = " << time << std::endl;
    std::cout << "======================================================================\n";
  }
  
}


