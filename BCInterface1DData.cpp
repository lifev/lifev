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
 *  @brief File containing the BCInterface1D_Data class
 *
 *  @date 10-05-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */


#include <lifemc/lifesolver/BCInterface1DData.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterface1D_Data::BCInterface1D_Data() :
        M_side               (),
        M_line               (),
        M_quantity           (),
        M_baseString         (),
        M_base               (),
        M_resistance         (),
        M_mapSide            (),
        M_mapQuantity        (),
        M_mapLine            (),
        M_mapBase            ()
{
    //Set mapSide
    M_mapSide["left"]   = OneD_left;
    M_mapSide["right"]  = OneD_right;

    //Set mapQuantity
    M_mapQuantity["A"]      = OneD_A;
    M_mapQuantity["Q"]      = OneD_Q;
    M_mapQuantity["W1"]     = OneD_W1;
    M_mapQuantity["W2"]     = OneD_W2;
    M_mapQuantity["P"]      = OneD_P;

    //Set mapLine
    M_mapLine["first"]  = OneD_first;
    M_mapLine["second"] = OneD_second;

    //Set mapBase
    M_mapBase["function"]         = BCInterface1D_function;
    M_mapBase["functionFile"]     = BCInterface1D_functionFile;
    M_mapBase["OPERfunction"]     = BCInterface1D_OPERfunction;
    M_mapBase["OPERfunctionFile"] = BCInterface1D_OPERfunctionFile;
    M_mapBase["Default"]          = BCInterface1D_Default;
}

BCInterface1D_Data::BCInterface1D_Data( const BCInterface1D_Data& data ) :
        M_side              ( data.M_side ),
        M_line              ( data.M_line ),
        M_quantity          ( data.M_quantity ),
        M_baseString        ( data.M_baseString ),
        M_base              ( data.M_base ),
        M_resistance        ( data.M_resistance ),
        M_mapSide           ( data.M_mapSide ),
        M_mapQuantity       ( data.M_mapQuantity ),
        M_mapLine           ( data.M_mapLine ),
        M_mapBase           ( data.M_mapBase )
{
}

// ===================================================
// Operators
// ===================================================
BCInterface1D_Data&
BCInterface1D_Data::operator=( const BCInterface1D_Data& data )
{
    if ( this != &data )
    {
        M_side              = data.M_side;
        M_line              = data.M_line;
        M_quantity          = data.M_quantity;
        M_baseString        = data.M_baseString;
        M_base              = data.M_base;
        M_resistance        = data.M_resistance;
        M_mapSide           = data.M_mapSide;
        M_mapQuantity       = data.M_mapQuantity;
        M_mapLine           = data.M_mapLine;
        M_mapBase           = data.M_mapBase;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
BCInterface1D_Data::readBC( const std::string& fileName, const std::string& dataSection, const BCName& name )
{
    GetPot dataFile( fileName );

    readSide( dataFile, ( dataSection + name + "/side" ).c_str() );
    readQuantity( dataFile, ( dataSection + name + "/quantity" ).c_str() );
    readLine( dataFile, ( dataSection + name + "/line" ).c_str() );
    readBase( dataFile, dataSection + name + "/" );
    readResistance( dataFile, ( dataSection + name + "/resistance" ).c_str() );
}

void
BCInterface1D_Data::showMe( std::ostream& output ) const
{
    output << "Side       = " << M_side << std::endl;
    output << "Line       = " << M_line << std::endl;
    output << "Quantity   = " << M_quantity << std::endl;
    output << "base       = " << M_base.second << std::endl;
    output << "baseString = " << M_baseString << std::endl;
    output << "Resistance  = ";
    for ( UInt i(0); i < static_cast<UInt>( M_resistance.size() ); ++i )
        output << M_resistance[i] << " ";
    output << "\n";
}

// ===================================================
// Methods
// ===================================================
void
BCInterface1D_Data::setBaseString( const std::string& baseString )
{
    M_baseString = baseString;
    boost::replace_all( M_baseString, " ", "" );
}

// ===================================================
// Private Methods
// ===================================================
void
BCInterface1D_Data::readBase( const GetPot& dataFile, const std::string& base )
{
    for ( std::map< std::string, bcBaseList_Type >::iterator j = M_mapBase.begin(); j
            != M_mapBase.end(); ++j )
        if ( isBase( dataFile, ( base + j->first ).c_str() ) )
        {
            M_base.first = j->first;
            M_base.second = M_mapBase[j->first];

            break;
        }
}

bool
BCInterface1D_Data::isBase( const GetPot& dataFile, const char* base )
{
    M_baseString = dataFile( base, " " );

    return dataFile.checkVariable( base );
}

void
BCInterface1D_Data::readResistance( const GetPot& dataFile, const char* resistance )
{
    UInt resistanceSize = dataFile.vector_variable_size( resistance );

    M_resistance.clear();
    M_resistance.reserve( resistanceSize );

    for ( UInt j( 0 ); j < resistanceSize; ++j )
        M_resistance.push_back( dataFile( resistance, 0, j ) );
}

} // Namespace LifeV
