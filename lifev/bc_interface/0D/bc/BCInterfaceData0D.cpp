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
 *  @file
 *  @brief File containing the BCInterfaceData0D class
 *
 *  @date 17-07-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/bc_interface/fem/BCInterfaceData0D.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterfaceData0D::BCInterfaceData0D() :
    BCInterfaceData         (),
    M_flag                  (),
    M_type                  (),
    M_mapType               ()
{
    //Set mapType
    M_mapType["Current"]  = Current;
    M_mapType["Voltage"]  = Voltage;
}

BCInterfaceData0D::BCInterfaceData0D ( const BCInterfaceData0D& data ) :
    BCInterfaceData         ( data ),
    M_flag                  ( data.M_flag ),
    M_type                  ( data.M_type ),
    M_mapType               ( data.M_mapType )
{
}

// ===================================================
// Operators
// ===================================================
BCInterfaceData0D&
BCInterfaceData0D::operator= ( const BCInterfaceData0D& data )
{
    if ( this != &data )
    {
        BCInterfaceData::operator= ( data );
        M_flag                  = data.M_flag;
        M_type                  = data.M_type;
        M_mapType               = data.M_mapType;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
BCInterfaceData0D::readBC ( const std::string& fileName, const std::string& dataSection, const std::string& name )
{
    // Call to the base class
    dataContainer_Type::readBC ( fileName, dataSection, name );

    // Read 0D data
    GetPot dataFile ( fileName );

    readFlag ( dataFile, ( dataSection + name + "/flag" ).c_str() );
    readType ( dataFile, ( dataSection + name + "/type0D" ).c_str() );
}

void
BCInterfaceData0D::showMe ( std::ostream& output ) const
{
    // Call to the base class
    dataContainer_Type::showMe ( output );

    // Show 0D data
    output << "Flag              = " << static_cast< Real > ( M_flag ) << std::endl;
    output << "Type              = " << M_type << std::endl;
}

} // Namespace LifeV
