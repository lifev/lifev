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
 *  @brief File containing the BCInterfaceData class
 *
 *  @date 17-07-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <life/lifesolver/BCInterfaceData.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterfaceData::BCInterfaceData() :
        M_base                  (),
        M_baseString            (),
        M_baseRobinAlpha        (),
        M_baseStringRobinAlpha  (),
        M_baseRobinBeta         (),
        M_baseStringRobinBeta   (),
        M_baseDirectional       (),
        M_baseStringDirectional (),
        M_mapBase               (),
        M_side                  (),
        M_line                  (),
        M_quantity              (),
        M_resistance            (),
        M_capacitance           (),
        M_mapSide               (),
        M_mapQuantity           (),
        M_mapLine               (),
        M_name                  (),
        M_flag                  (),
        M_type                  (),
        M_mode                  (),
        M_componentsVector      (),
        M_mapType               (),
        M_mapMode               ()
{
    //Set mapBase
    M_mapBase["function"]           = BCIFunction;
    M_mapBase["functionFile"]       = BCIFunctionFile;
    M_mapBase["functionSolver"]     = BCIFunctionSolver;
    M_mapBase["functionFileSolver"] = BCIFunctionFileSolver;
    M_mapBase["functionDefault"]    = BCI1DFunctionDefault;
    M_mapBase["dataInterpolator"]   = BCI3DDataInterpolator;
    M_mapBase["FSI"]                = BCI3DFSI;

    //Set mapSide
    M_mapSide["left"]   = OneDimensional::left;
    M_mapSide["right"]  = OneDimensional::right;

    //Set mapQuantity
    M_mapQuantity["A"]  = OneDimensional::A;
    M_mapQuantity["Q"]  = OneDimensional::Q;
    M_mapQuantity["W1"] = OneDimensional::W1;
    M_mapQuantity["W2"] = OneDimensional::W2;
    M_mapQuantity["P"]  = OneDimensional::P;
    M_mapQuantity["S"]  = OneDimensional::S;

    //Set mapLine
    M_mapLine["first"]  = OneDimensional::first;
    M_mapLine["second"] = OneDimensional::second;

    //Set mapType
    M_mapType["Essential"]         = Essential;
    M_mapType["EssentialEdges"]    = EssentialEdges;
    M_mapType["EssentialVertices"] = EssentialVertices;
    M_mapType["Natural"]           = Natural;
    M_mapType["Robin"]             = Robin;
    M_mapType["Flux"]              = Flux;
    M_mapType["Resistance"]        = Resistance;

    //Set mapMode
    M_mapMode["Scalar"]      = Scalar;
    M_mapMode["Full"]        = Full;
    M_mapMode["Component"]   = Component;
    M_mapMode["Normal"]      = Normal;
    M_mapMode["Tangential"]  = Tangential;
    M_mapMode["Directional"] = Directional;
}

BCInterfaceData::BCInterfaceData( const BCInterfaceData& data ) :
        M_base                  ( data.M_base ),
        M_baseString            ( data.M_baseString ),
        M_baseRobinAlpha        ( data.M_baseRobinAlpha ),
        M_baseStringRobinAlpha  ( data.M_baseStringRobinAlpha ),
        M_baseRobinBeta         ( data.M_baseRobinBeta ),
        M_baseStringRobinBeta   ( data.M_baseStringRobinBeta ),
        M_baseDirectional       ( data.M_baseDirectional ),
        M_baseStringDirectional ( data.M_baseStringDirectional ),
        M_mapBase               ( data.M_mapBase ),
        M_side                  ( data.M_side ),
        M_line                  ( data.M_line ),
        M_quantity              ( data.M_quantity ),
        M_resistance            ( data.M_resistance ),
        M_capacitance           ( data.M_capacitance ),
        M_mapSide               ( data.M_mapSide ),
        M_mapQuantity           ( data.M_mapQuantity ),
        M_mapLine               ( data.M_mapLine ),
        M_name                  ( data.M_name ),
        M_flag                  ( data.M_flag ),
        M_type                  ( data.M_type ),
        M_mode                  ( data.M_mode ),
        M_componentsVector      ( data.M_componentsVector ),
        M_mapType               ( data.M_mapType ),
        M_mapMode               ( data.M_mapMode )
{
}

// ===================================================
// Operators
// ===================================================
BCInterfaceData&
BCInterfaceData::operator=( const BCInterfaceData& data )
{
    if ( this != &data )
    {
        M_base                  = data.M_base;
        M_baseString            = data.M_baseString;
        M_baseRobinAlpha        = data.M_baseRobinAlpha;
        M_baseStringRobinAlpha  = data.M_baseStringRobinAlpha;
        M_baseRobinBeta         = data.M_baseRobinBeta;
        M_baseStringRobinBeta   = data.M_baseStringRobinBeta;
        M_baseDirectional       = data.M_baseDirectional;
        M_baseStringDirectional = data.M_baseStringDirectional;
        M_mapBase               = data.M_mapBase;
        M_side                  = data.M_side;
        M_line                  = data.M_line;
        M_quantity              = data.M_quantity;
        M_resistance            = data.M_resistance;
        M_capacitance           = data.M_capacitance;
        M_mapSide               = data.M_mapSide;
        M_mapQuantity           = data.M_mapQuantity;
        M_mapLine               = data.M_mapLine;
        M_name                  = data.M_name;
        M_flag                  = data.M_flag;
        M_type                  = data.M_type;
        M_mode                  = data.M_mode;
        M_componentsVector      = data.M_componentsVector;
        M_mapType               = data.M_mapType;
        M_mapMode               = data.M_mapMode;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
BCInterfaceData::readBC( const std::string& fileName, const std::string& dataSection, const bcName_Type& name )
{
    GetPot dataFile( fileName );

    // 1D
    readSide( dataFile, ( dataSection + name + "/side" ).c_str() );
    readQuantity( dataFile, ( dataSection + name + "/quantity" ).c_str() );
    readLine( dataFile, ( dataSection + name + "/line" ).c_str() );
    readResistance( dataFile, ( dataSection + name + "/resistance" ).c_str() );
    readCapacitance( dataFile, ( dataSection + name + "/capacitance" ).c_str() );

    // 3D
    M_name = name;

    readFlag( dataFile, ( dataSection + name + "/flag" ).c_str() );
    readType( dataFile, ( dataSection + name + "/type" ).c_str() );
    readMode( dataFile, ( dataSection + name + "/mode" ).c_str() );
    readComV( dataFile, ( dataSection + name + "/component" ).c_str() );

    // Read base
    readBase( dataFile, dataSection + name + "/", M_base, M_baseString );
    if ( M_type == Robin )
    {
        readBase( dataFile, dataSection + name + "/RobinAlpha/", M_baseRobinAlpha, M_baseStringRobinAlpha );
        readBase( dataFile, dataSection + name + "/RobinBeta/", M_baseRobinBeta, M_baseStringRobinBeta );
    }
    if ( M_mode == Directional )
        readBase( dataFile, dataSection + name + "/Directional/", M_baseDirectional, M_baseStringDirectional );
}

void
BCInterfaceData::showMe( std::ostream& output ) const
{
    output << "baseString = " << M_baseString << std::endl;

    output << "Side       = " << M_side << std::endl;
    output << "Line       = " << M_line << std::endl;
    output << "Quantity   = " << M_quantity << std::endl;
    output << "Resistance  = ";
    for ( UInt i(0); i < static_cast<UInt>( M_resistance.size() ); ++i )
        output << M_resistance[i] << " ";
    output << "\n";
    output << "Capacitance= " << M_capacitance << std::endl;
    output << "base       = " << M_base.second << std::endl;

    output << "Flag       = " << static_cast< Real > ( M_flag ) << std::endl;
    output << "Type       = " << M_type << std::endl;
    output << "Mode       = " << M_mode << std::endl;
    output << "comV:      = ";
    for ( UInt i(0); i < static_cast<UInt>( M_componentsVector.size() ); ++i )
        output << M_componentsVector[i] << " ";
    output << "\n";
}

// ===================================================
// Set Methods
// ===================================================
void BCInterfaceData::setBaseString( const std::string& baseString )
{
    M_baseString = baseString;
    boost::replace_all( M_baseString, " ", "" );
}

// ===================================================
// Private Methods
// ===================================================
void
BCInterfaceData::readResistance( const GetPot& dataFile, const char* resistance )
{
    UInt resistanceSize = dataFile.vector_variable_size( resistance );

    M_resistance.resize( resistanceSize );
    for ( UInt j( 0 ); j < resistanceSize; ++j )
        M_resistance[j] = dataFile( resistance, 0, j );
}

void
BCInterfaceData::readComV( const GetPot& dataFile, const char* component )
{
    UInt componentSize = dataFile.vector_variable_size( component );

    M_componentsVector.resize( componentSize );
    for ( UInt j( 0 ); j < componentSize; ++j )
        M_componentsVector[j] = dataFile( component, 0, j );
}

void
BCInterfaceData::readBase( const GetPot& dataFile, const std::string& path, std::pair< std::string, baseList_Type >& base, std::string& baseString )
{
    for ( std::map< std::string, baseList_Type >::iterator j = M_mapBase.begin(); j != M_mapBase.end(); ++j )
        if ( isBase( dataFile, ( path + j->first ).c_str(), baseString ) )
        {
            base.first = j->first;
            base.second = M_mapBase[j->first];

            break;
        }
}

bool
BCInterfaceData::isBase( const GetPot& dataFile, const char* base, std::string& baseString )
{
    baseString = dataFile( base, " " );

    return dataFile.checkVariable( base );
}

} // Namespace LifeV
