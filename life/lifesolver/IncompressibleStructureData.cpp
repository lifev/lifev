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
    @brief File containing the implementation of the file IncompressibleStructureData.hpp

    @author Simone Rossi <simone.rossi@epfl.ch>
    @contributor
    @maintainer Simone Rossi <simone.rossi@epfl.ch>

    @date 04-2011

 */


#include <life/lifesolver/IncompressibleStructureData.hpp>
#include <life/lifecore/LifeV.hpp>


namespace LifeV
{


// ===================================================
// Constructors
// ===================================================

IncompressibleStructureData::IncompressibleStructureData() :
        M_viscosity                        ( ),
        M_uOrder                           ( ),
        M_pOrder                           ( ),
        M_verbose                          ( )
{
}



IncompressibleStructureData::IncompressibleStructureData( const IncompressibleStructureData& incompressibleStructureData ) :
        M_viscosity                        ( incompressibleStructureData.M_viscosity ),
        M_uOrder                           ( incompressibleStructureData.M_uOrder ),
        M_pOrder                           ( incompressibleStructureData.M_pOrder ),
        M_verbose                          ( incompressibleStructureData.M_verbose )
{
}






// ===================================================
// Methods
// ===================================================

IncompressibleStructureData&
IncompressibleStructureData::operator=( const IncompressibleStructureData& incompressibleStructureData )
{
    if ( this != &incompressibleStructureData )
    {
        M_viscosity                        = incompressibleStructureData.M_viscosity;
        M_uOrder                           = incompressibleStructureData.M_uOrder;
        M_pOrder                           = incompressibleStructureData.M_pOrder;
        M_verbose                          = incompressibleStructureData.M_verbose;
    }

    return *this;
}


void
IncompressibleStructureData::setup( const GetPot& dataFile, const std::string& section )
{
    // Physics
    M_viscosity.push_back ( dataFile( ( section + "/physics/viscosity" ).data(), 1. ) );

    // FE Order
    M_uOrder       = dataFile( ( section + "/space_discretization/disp_order" ).data(), "P1");
    M_pOrder       = dataFile( ( section + "/space_discretization/press_order" ).data(), "P1");

    // Miscellaneous
    M_verbose      = 1; //dataFile( ( section + "/miscellaneous/verbose" ).data(), 1 );

}


void
IncompressibleStructureData::showMe( std::ostream& output ) const
{
    output << "\n*** Values for data [structure/physics]\n\n";
    output << "viscosity   = " << M_viscosity[0] << std::endl;
    output << std::endl;


}


} //end namespace LifeV
