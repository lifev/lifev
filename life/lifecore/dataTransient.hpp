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
/*!
  \file dataTransient.hpp
  \author J.F. Gerbeau
  \date 09/2004
  \version 1.0
 
  \brief File containing a class for handling temporal discretization with GetPot
 
*/
#ifndef _DATATRANSIENT_H_
#define _DATATRANSIENT_H_
#include <string>
#include <iostream>
#include "GetPot.hpp"


namespace LifeV
{

/*
  \author JFG
  \brief very poor data for time dependent problem
 
  \todo merge with dataTime
 
  \todo allow variable time steps
  \todo select a stopping test (based on either max_time_iter or max_time)
  \todo a standard banner for new time step
  \todo tolerance for steady state
 
*/ 
//using namespace std;
class DataTransient
{
public:
    double max_time;
    int max_time_iter;
    double time_step;
    int post_proc_period;
    /**
       Constructor from GetPot
     */
    DataTransient( const GetPot& dfile );
    /**
       Print information
     */
    void dataTransientShowMe( std::ostream& c );
    /**
       Print some help
     */
    void dataTransientHelp( std::ostream& c );

    /**
       Print current iteration and time
     */
    void timeBanner( int iter, double t );
};
}
#endif
