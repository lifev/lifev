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
  \file dataNavierStokes.h
  \author M.A. Fernandez
  \date 01/2003 
  \version 1.0

  \brief File containing a class for handling NavierStokes data with GetPot

*/
#ifndef _DATANAVIERSTOKES_H_
#define _DATANAVIERSTOKES_H_
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "lifeV.hpp"
#include "dataMesh.hpp"
#include "dataTime.hpp"

/*! 
  \class DataNavierStokes

  Base class which holds usual data for the NavierStokes equations solvers

*/
template <typename Mesh>
class DataNavierStokes:
public DataMesh<Mesh>,
public DataTime {
 public:

  //! Constructor
  DataNavierStokes(const GetPot& dfile);
  
  //! Ouptut
  void showMe(ostream& c=cout);
  //! End time
  Real density() const;
  Real viscosity() const;
  Real inittime() const;
  Real endtime() const;

  UInt verbose() const;
  Real dump_init() const;
  UInt dump_period() const;

 protected:
  //! Physics
  Real _rho; // density
  Real _mu; // viscosity
  Real _inittime; // initial time (Alex December 2003)
  Real _endtime; // end time


  //! Miscellaneous 
  UInt _verbose; // temporal output verbose
  Real _dump_init; // time for starting the dumping of the results (Alex December 2003)
  UInt _dump_period; // frequency of the dumping (one dump after _dump_period time steps) (Alex December 2003)
};


//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh>
DataNavierStokes<Mesh>::
DataNavierStokes(const GetPot& dfile):
  DataMesh<Mesh>(dfile,"fluid/discretization"),
  DataTime(dfile,"fluid/discretization") {
  
  // physics
  _rho       = dfile("fluid/physics/density",1.);
  _mu        = dfile("fluid/physics/viscosity",1.);
  _inittime  = dfile("fluid/physics/inittime",0.);
  _endtime   = dfile("fluid/physics/endtime",1.);

  //miscellaneous
  _verbose   = dfile("fluid/miscellaneous/verbose",1);
  _dump_init = dfile("fluid/miscellaneous/dump_init",_inittime);
  _dump_period = dfile("fluid/miscellaneous/dump_period",1);
}

// Output
template <typename Mesh>
void DataNavierStokes<Mesh>::
showMe(ostream& c)
{
  // physics
  c << "\n*** Values for data [fluid/physics]\n\n";
  c << "density   = " << _rho << endl; 
  c << "viscosity = " << _mu << endl;
  c << "initial time = " << _inittime << endl;
  c << "endtime   = " << _endtime << endl;

  c << "\n*** Values for data [fluid/miscellaneous]\n\n";
  c << "verbose   = " << _verbose << endl;
  c << "initial time for writing solution  = " << _dump_init << endl;
  c << "number of time steps between two consecutive dumps of the solution = " << _dump_period << endl;

  c << "\n*** Values for data [fluid/discretization]\n\n";
  DataMesh<Mesh>::showMe(c);
  DataTime::showMe(c);

}
////////////////////
// The density
template <typename Mesh>
Real DataNavierStokes<Mesh>::
density() const {
  return  _density;
}

// The viscosity
template <typename Mesh>
Real DataNavierStokes<Mesh>::
viscosity() const {
  return  _viscosity;
}


// The initial time
template <typename Mesh>
Real DataNavierStokes<Mesh>::
inittime() const {
  return  _inittime;
}

// The end time
template <typename Mesh>
Real DataNavierStokes<Mesh>::
endtime() const {
  return  _endtime;
}
////////////////

// verbose variable
template <typename Mesh>
UInt DataNavierStokes<Mesh>::
verbose() const {
  return  _verbose;
}
 // Dumping start
template <typename Mesh>
Real DataNavierStokes<Mesh>::
dump_init() const {
  return  _dump_init;
}
// Period of dumping
template <typename Mesh>
UInt DataNavierStokes<Mesh>::
dump_period() const {
  return  _dump_period;
}

#endif
