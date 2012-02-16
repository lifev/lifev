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
    @brief File containing a class for handling temporal discretization

    @author M.A. Fernandez
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @date 01-06-2009

    @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>

    @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 */
#include <life/lifefem/TimeData.hpp>

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
        M_timeStepNumber( 0 ),
        M_orderBDF      ( 1 ),
        M_theta         ( 0.25 ),
        M_gamma         ( 0.5)
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
    M_periodTime    = timeData.M_periodTime;
    M_time            = timeData.M_time;
    M_timeStep        = timeData.M_timeStep;
    M_timeStepNumber= timeData.M_timeStepNumber;
    M_orderBDF        = timeData.M_orderBDF;
    M_theta            = timeData.M_theta;
    M_gamma         = timeData.M_gamma;
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
    M_orderBDF = dataFile(( section + "/BDF_order" ).data(), 1 );
    M_theta = dataFile((section + "/theta").data(),0.25);
    M_gamma = dataFile(( section + "/zeta").data(),0.5);
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
    output << "BDF order      = " << M_orderBDF       << std::endl;
    output << "theta          = " << M_theta          << std::endl;
    output << "gamma          = " << M_gamma          << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
std::vector<Real>
TimeData::coefficientsNewmark()
{
    std::vector<Real> coefficients;

    coefficients.push_back( M_theta );
    coefficients.push_back( M_gamma );

    return  coefficients;
}

// ===================================================
// Private Methods
// ===================================================
Real
TimeData::round( const Real& n, const Int& decimal ) const
{
    return std::floor( n * std::pow(10.0, decimal) + 0.5 )  / std::pow(10.0, decimal);
}

}
