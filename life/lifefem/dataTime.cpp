/*
This file is part of the LifeV library
Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
  \file dataTime.cpp

  \version 1.0
  \date 01/2003
  \author M.A. Fernandez

  \version 1.9
  \date 06/2009
  \author Cristiano Malossi<cristiano.malossi@epfl.ch>

  \brief File containing a class for handling temporal discretization
*/
#include <life/lifefem/dataTime.hpp>

namespace LifeV
{

// Constructors
DataTime::DataTime( const GetPot& dfile, const std::string& section ) :
	M_initialTime	( dfile(( section + "/initialtime" 	).data(), 0.) ),
	M_endTime		( dfile(( section + "/endtime" 		).data(), 1.) ),
	M_time			( 0. ),
	M_timeStep		( dfile(( section + "/timestep" ).data(), 1.) ),
	M_BDF_order		( dfile(( section + "/BDF_order" ).data(), 1 ) )
{
}

DataTime::DataTime(const DataTime& dataTime)
{
	M_initialTime	= dataTime.M_initialTime;
	M_endTime		= dataTime.M_endTime;
	M_time			= dataTime.M_time;
	M_timeStep		= dataTime.M_timeStep;
	M_BDF_order		= dataTime.M_BDF_order;
}

// Output
void DataTime::showMe( std::ostream& output ) const
{
	output << "Initial time = " << M_initialTime	<< std::endl;
	output << "End time     = " << M_endTime		<< std::endl;
	output << "Time         = " << M_time			<< std::endl;
	output << "Timestep     = " << M_timeStep		<< std::endl;
	output << "BDF order    = " << M_BDF_order  	<< std::endl;
}

}
