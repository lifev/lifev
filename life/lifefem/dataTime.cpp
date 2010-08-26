//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a class for handling temporal discretization
 *
 *  @author M.A. Fernandez
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 01-06-2009
 */
#include <life/lifefem/dataTime.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
DataTime::DataTime( ) :
    M_initialTime   ( 0. ),
    M_endTime       ( 1. ),
    M_time          ( M_initialTime ),
    M_timeStep      ( M_endTime ),
    M_BDF_order     ( 1 ),
    M_theta         ( 0.25 ),
    M_zeta          ( 0.5)
{
}

DataTime::DataTime( const GetPot& dfile, const std::string& section )
{
    setup( dfile, section );
}

DataTime::DataTime( const DataTime& dataTime )
{
	M_initialTime	= dataTime.M_initialTime;
	M_endTime		= dataTime.M_endTime;
	M_time			= dataTime.M_time;
	M_timeStep		= dataTime.M_timeStep;
	M_BDF_order		= dataTime.M_BDF_order;
	M_theta	        = dataTime.M_theta;
	M_zeta          = dataTime.M_zeta;
}

// ===================================================
// Methods
// ===================================================
void
DataTime::setup( const GetPot& dfile, const std::string& section )
{
    M_initialTime = dfile(( section + "/initialtime"  ).data(), 0.);
    M_endTime = dfile(( section + "/endtime"      ).data(), 1.);
    M_time = M_initialTime;
    M_timeStep = dfile(( section + "/timestep" ).data(), M_endTime );
    M_BDF_order = dfile(( section + "/BDF_order" ).data(), 1 );
    M_theta = dfile((section + "/theta").data(),0.25);
    M_zeta = dfile(( section + "/zeta").data(),0.5);
}

void
DataTime::showMe( std::ostream& output ) const
{
    output << "\n*** DataTime: values for user-defined data\n";

    output << "\n[/time_discretization]" << std::endl;
    output << "Initial time = " << M_initialTime << std::endl;
	output << "End time     = " << M_endTime	 << std::endl;
	output << "Time         = " << M_time		 << std::endl;
	output << "TimeStep     = " << M_timeStep	 << std::endl;
	output << "BDF order    = " << M_BDF_order   << std::endl;
	output << "theta        = " << M_theta       << std::endl;
	output << "zeta         = " << M_zeta        << std::endl;
}

// ===================================================
// Methods
// ===================================================
std::vector<Real>
DataTime::getNewmark_parameters()
{
    std::vector<Real> parameters;

    parameters.push_back( M_theta );
    parameters.push_back( M_zeta );

    return parameters;
}

// ===================================================
// Private Methods
// ===================================================
Real
DataTime::round( const Real n, const Int decimal ) const
{
    return std::floor( n * std::pow(10.0, decimal) + 0.5 )  / std::pow(10.0, decimal);
}

}
