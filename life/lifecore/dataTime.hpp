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


namespace LifeV
{

/*!
  \class DataTime

  Base class which holds data concerning temporal discretization

*/
class DataTime
{
public:

    //! Constructor
    DataTime( const GetPot& dfile, const std::string& section = "discretization" );

    //! Ouptut
    virtual void showMe( std::ostream& c = std::cout ) const;

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
}
#endif
