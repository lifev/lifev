/*!
  \file dataTime.h
  \author M.A. Fernandez
  \date 01/2003 
  \version 1.0

  \brief File containing a class for handling temporal discretization with GetPot

*/
#ifndef _DATATIME_H_
#define _DATATIME_H_
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "lifeV.hpp"


/*! 
  \class DataMesh

  Base class which holds data concerning temporal discretization

*/
class DataTime
{
 public:

  //! Constructor
  DataTime(const GetPot& dfile, const string& section="discretization");
  
  //! Ouptut
  virtual void showMe(ostream& c=cout) const;

  //! Time step
  Real timestep() const;

  //! Order BDF formula
  unsigned int order_bdf() const;
  
  //! Virtual destructor
  virtual ~DataTime();

 protected:
 
  Real _dt; // time step
  unsigned int _order_bdf; //order of the time discretization formula
};

#endif
