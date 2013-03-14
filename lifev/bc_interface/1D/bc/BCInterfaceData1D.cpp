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
 *  @brief File containing the BCInterfaceData1D class
 *
 *  @date 17-07-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/bc_interface/1D/bc/BCInterfaceData1D.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterfaceData1D::BCInterfaceData1D() :
    BCInterfaceData         (),
    M_line                  (),
    M_quantity              (),
    M_resistance            (),
    M_capacitance           (),
    M_mapSide               (),
    M_mapQuantity           (),
    M_mapLine               ()
{
    //Set mapSide
    M_mapSide["left"]   = OneDFSI::left;
    M_mapSide["right"]  = OneDFSI::right;

    //Set mapQuantity
    M_mapQuantity["A"]  = OneDFSI::A;
    M_mapQuantity["Q"]  = OneDFSI::Q;
    M_mapQuantity["W1"] = OneDFSI::W1;
    M_mapQuantity["W2"] = OneDFSI::W2;
    M_mapQuantity["P"]  = OneDFSI::P;
    M_mapQuantity["S"]  = OneDFSI::S;

    //Set mapLine
    M_mapLine["first"]  = OneDFSI::first;
    M_mapLine["second"] = OneDFSI::second;
}

BCInterfaceData1D::BCInterfaceData1D ( const BCInterfaceData1D& data ) :
    BCInterfaceData         ( data ),
    M_line                  ( data.M_line ),
    M_quantity              ( data.M_quantity ),
    M_resistance            ( data.M_resistance ),
    M_capacitance           ( data.M_capacitance ),
    M_mapSide               ( data.M_mapSide ),
    M_mapQuantity           ( data.M_mapQuantity ),
    M_mapLine               ( data.M_mapLine )
{
}

// ===================================================
// Operators
// ===================================================
BCInterfaceData1D&
BCInterfaceData1D::operator= ( const BCInterfaceData1D& data )
{
    if ( this != &data )
    {
        BCInterfaceData::operator= ( data );
        M_line                  = data.M_line;
        M_quantity              = data.M_quantity;
        M_resistance            = data.M_resistance;
        M_capacitance           = data.M_capacitance;
        M_mapSide               = data.M_mapSide;
        M_mapQuantity           = data.M_mapQuantity;
        M_mapLine               = data.M_mapLine;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
BCInterfaceData1D::readBC ( const std::string& fileName, const std::string& dataSection, const std::string& name )
{
    // Call to the base class
    dataContainer_Type::readBC ( fileName, dataSection, name );

    // Read 3D data
    GetPot dataFile ( fileName );

    readSide ( dataFile, ( dataSection + name + "/side" ).c_str() );
    readQuantity ( dataFile, ( dataSection + name + "/quantity" ).c_str() );
    readLine ( dataFile, ( dataSection + name + "/line" ).c_str() );
    readResistance ( dataFile, ( dataSection + name + "/resistance" ).c_str() );
    readCapacitance ( dataFile, ( dataSection + name + "/capacitance" ).c_str() );
}

void
BCInterfaceData1D::showMe ( std::ostream& output ) const
{
    // Call to the base class
    dataContainer_Type::showMe ( output );

    // Show 1D data
    output << "Line              = " << M_line << std::endl;
    output << "Quantity          = " << M_quantity << std::endl;
    output << "Resistance        = ";
    for ( UInt i (0); i < static_cast<UInt> ( M_resistance.size() ); ++i )
    {
        output << M_resistance[i] << " ";
    }
    output << "\n";
    output << "Capacitance       = " << M_capacitance << std::endl;
}

// ===================================================
// Private Methods
// ===================================================
void
BCInterfaceData1D::readResistance ( const GetPot& dataFile, const char* resistance )
{
    UInt resistanceSize = dataFile.vector_variable_size ( resistance );

    M_resistance.resize ( resistanceSize );
    for ( UInt j ( 0 ); j < resistanceSize; ++j )
    {
        M_resistance[j] = dataFile ( resistance, 0, j );
    }
}

} // Namespace LifeV
