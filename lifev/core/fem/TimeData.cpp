//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief File containing a class for handling temporal data.

    @author M.A. Fernandez
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @date 01-06-2009

    @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 */
#include <lifev/core/fem/TimeData.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
TimeData::TimeData( ) :
        M_initialTime   ( 0. ),
        M_endTime       ( 1. ),
        M_periodTime    ( 1. ),
        M_time          ( M_initialTime ),
        M_timeStep      ( M_endTime ),
        M_timeStepNumber( 0 )
{
}

TimeData::TimeData( const GetPot& dataFile, const std::string& section ) :
        M_timeStepNumber( 0 )
{
    setup( dataFile, section );
}

TimeData::TimeData( const TimeData& timeData )
{
    M_initialTime    = timeData.M_initialTime;
    M_endTime        = timeData.M_endTime;
    M_periodTime     = timeData.M_periodTime;
    M_time           = timeData.M_time;
    M_timeStep       = timeData.M_timeStep;
    M_timeStepNumber = timeData.M_timeStepNumber;
}

// ===================================================
// Methods
// ===================================================
void
TimeData::setup( const GetPot& dataFile, const std::string& section )
{
    M_initialTime = dataFile(( section + "/initialtime"  ).data(), 0.);
    M_endTime = dataFile(( section + "/endtime"      ).data(), 1.);
    M_periodTime = dataFile(( section + "/periodtime"      ).data(), 1.);
    M_time = M_initialTime;
    M_timeStep = dataFile(( section + "/timestep" ).data(), M_endTime );
}

void
TimeData::showMe( std::ostream& output ) const
{
    output << "\n*** TimeData: values for user-defined data\n";

    output << "\n[/time_discretization]" << std::endl;
    output << "Initial time   = " << M_initialTime    << std::endl;
    output << "End time       = " << M_endTime          << std::endl;
    output << "Time           = " << M_time              << std::endl;
    output << "TimeStep       = " << M_timeStep          << std::endl;
    output << "TimeStepNumber = " << M_timeStepNumber << std::endl;
}



// ===================================================
// Utility methods
// ===================================================
Real
TimeData::round( const Real& n, const Int& decimal ) const
{
    return std::floor( n * std::pow(10.0, decimal) + 0.5 )  / std::pow(10.0, decimal);
}

}
