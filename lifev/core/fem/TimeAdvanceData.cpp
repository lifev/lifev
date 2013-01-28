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
    @brief File containing a class for handling temporal discretization.

    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @date 11-06-2012

    @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 */
#include <lifev/core/fem/TimeAdvanceData.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
TimeAdvanceData::TimeAdvanceData( ) :
        M_orderBDF      ( 1 ),
        M_theta         ( 0.25 ),
        M_gamma         ( 0.5 )
{
}

TimeAdvanceData::TimeAdvanceData( const GetPot& dataFile, const std::string& section )
{
    setup( dataFile, section );
}

TimeAdvanceData::TimeAdvanceData( const TimeAdvanceData& timeData )
{
    M_orderBDF       = timeData.M_orderBDF;
    M_theta          = timeData.M_theta;
    M_gamma          = timeData.M_gamma;
}

// ===================================================
// Methods
// ===================================================
void
TimeAdvanceData::setup( const GetPot& dataFile, const std::string& section )
{
    M_orderBDF = dataFile(( section + "/BDF_order" ).data(), 1 );
    M_theta = dataFile((section + "/theta").data(),0.25);
    M_gamma = dataFile(( section + "/gamma").data(),0.5);
}

void
TimeAdvanceData::showMe( std::ostream& output ) const
{
    output << "\n*** TimeAdvanceData: values for user-defined data\n";

    output << "\n[/time_discretization]" << std::endl;
    output << "BDF order      = " << M_orderBDF       << std::endl;
    output << "theta          = " << M_theta          << std::endl;
    output << "gamma          = " << M_gamma          << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
std::vector<Real>
TimeAdvanceData::coefficientsNewmark()
{
    std::vector<Real> coefficients;

    coefficients.push_back( M_theta );
    coefficients.push_back( M_gamma );

    return  coefficients;
}

}
