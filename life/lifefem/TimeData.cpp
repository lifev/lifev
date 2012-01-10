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

TimeData::TimeData( const GetPot& dfile, const std::string& section ) :
        M_timeStepNumber( 0 )
{
    setup( dfile, section );
}

TimeData::TimeData( const TimeData& TimeData )
{
    M_initialTime	= TimeData.M_initialTime;
    M_endTime		= TimeData.M_endTime;
    M_periodTime    = TimeData.M_periodTime;
    M_time			= TimeData.M_time;
    M_timeStep		= TimeData.M_timeStep;
    M_timeStepNumber= TimeData.M_timeStepNumber;
    M_orderBDF		= TimeData.M_orderBDF;
    M_theta	        = TimeData.M_theta;
    M_gamma         = TimeData.M_gamma;
}

// ===================================================
// Methods
// ===================================================
void
TimeData::setup( const GetPot& dfile, const std::string& section )
{
    M_initialTime = dfile(( section + "/initialtime"  ).data(), 0.);
    M_endTime = dfile(( section + "/endtime"      ).data(), 1.);
    M_periodTime = dfile(( section + "/periodtime"      ).data(), 1.);
    M_time = M_initialTime;
    M_timeStep = dfile(( section + "/timestep" ).data(), M_endTime );
    M_orderBDF = dfile(( section + "/BDF_order" ).data(), 1 );
    M_theta = dfile((section + "/theta").data(),0.25);
    M_gamma = dfile(( section + "/zeta").data(),0.5);
}

void
TimeData::showMe( std::ostream& output ) const
{
    output << "\n*** TimeData: values for user-defined data\n";

    output << "\n[/time_discretization]" << std::endl;
    output << "Initial time   = " << M_initialTime    << std::endl;
    output << "End time       = " << M_endTime	      << std::endl;
    output << "Time           = " << M_time		      << std::endl;
    output << "TimeStep       = " << M_timeStep	      << std::endl;
    output << "TimeStepNumber = " << M_timeStepNumber << std::endl;
    output << "BDF order      = " << M_orderBDF       << std::endl;
    output << "theta          = " << M_theta          << std::endl;
    output << "gamma          = " << M_gamma          << std::endl;
}

// ===================================================
// Methods
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
TimeData::round( const Real n, const Int decimal ) const
{
    return std::floor( n * std::pow(10.0, decimal) + 0.5 )  / std::pow(10.0, decimal);
}

}
