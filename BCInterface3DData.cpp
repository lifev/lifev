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
 *  @brief File containing the BCInterface3DData class
 *
 *  @date 17-07-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/BCInterface3DData.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterfaceData::BCInterfaceData() :
        M_baseString         (),
        M_side               (),
        M_line               (),
        M_quantity           (),
        M_base1D             (),
        M_resistance         (),
        M_mapSide            (),
        M_mapQuantity        (),
        M_mapLine            (),
        M_mapBase1D          (),
        M_name               (),
        M_flag               (),
        M_type               (),
        M_mode               (),
        M_comV               (),
        M_direction          (),
        M_base3D             (),
        M_mapType            (),
        M_mapMode            (),
        M_mapBase3D          ()
{
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

    //Set mapBase
    M_mapBase1D["function"]         = BCI1DFunction;
    M_mapBase1D["functionFile"]     = BCI1DFunctionFile;
    M_mapBase1D["OPERfunction"]     = BCI1DFunctionSolver;
    M_mapBase1D["OPERfunctionFile"] = BCI1DFunctionFileSolver;
    M_mapBase1D["Default"]          = BCI1DFunctionDefault;

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

    //Set mapBase
    M_mapBase3D["function"]         = BCI3DFunction;
    M_mapBase3D["functionFile"]     = BCI3DFunctionFile;
    M_mapBase3D["OPERfunction"]     = BCI3DFunctionSolver;
    M_mapBase3D["OPERfunctionFile"] = BCI3DFunctionFileSolver;
    M_mapBase3D["FSI"]              = BCI3DFunctionFSI;
}

BCInterfaceData::BCInterfaceData( const BCInterfaceData& data ) :
        M_baseString        ( data.M_baseString ),
        M_side              ( data.M_side ),
        M_line              ( data.M_line ),
        M_quantity          ( data.M_quantity ),
        M_base1D            ( data.M_base1D ),
        M_resistance        ( data.M_resistance ),
        M_mapSide           ( data.M_mapSide ),
        M_mapQuantity       ( data.M_mapQuantity ),
        M_mapLine           ( data.M_mapLine ),
        M_mapBase1D         ( data.M_mapBase1D ),
        M_name              ( data.M_name ),
        M_flag              ( data.M_flag ),
        M_type              ( data.M_type ),
        M_mode              ( data.M_mode ),
        M_comV              ( data.M_comV ),
        M_direction         ( data.M_direction ),
        M_base3D            ( data.M_base3D ),
        M_mapType           ( data.M_mapType ),
        M_mapMode           ( data.M_mapMode ),
        M_mapBase3D         ( data.M_mapBase3D )
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
        M_baseString        = data.M_baseString;
        M_side              = data.M_side;
        M_line              = data.M_line;
        M_quantity          = data.M_quantity;
        M_base1D            = data.M_base1D;
        M_resistance        = data.M_resistance;
        M_mapSide           = data.M_mapSide;
        M_mapQuantity       = data.M_mapQuantity;
        M_mapLine           = data.M_mapLine;
        M_mapBase1D         = data.M_mapBase1D;
        M_name              = data.M_name;
        M_flag              = data.M_flag;
        M_type              = data.M_type;
        M_mode              = data.M_mode;
        M_comV              = data.M_comV;
        M_direction         = data.M_direction;
        M_base3D            = data.M_base3D;
        M_mapType           = data.M_mapType;
        M_mapMode           = data.M_mapMode;
        M_mapBase3D         = data.M_mapBase3D;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
BCInterfaceData::readBC1D( const std::string& fileName, const std::string& dataSection, const bcName_Type& name )
{
    GetPot dataFile( fileName );

    readSide( dataFile, ( dataSection + name + "/side" ).c_str() );
    readQuantity( dataFile, ( dataSection + name + "/quantity" ).c_str() );
    readLine( dataFile, ( dataSection + name + "/line" ).c_str() );
    readBase1D( dataFile, dataSection + name + "/" );
    readResistance( dataFile, ( dataSection + name + "/resistance" ).c_str() );
}

void
BCInterfaceData::readBC3D( const std::string& fileName,
                         const std::string& dataSection,
                         const bcName_Type& name )
{
    GetPot dataFile( fileName );

    M_name = name;

    readFlag( dataFile, ( dataSection + name + "/flag" ).c_str() );
    readType( dataFile, ( dataSection + name + "/type" ).c_str() );
    readMode( dataFile, ( dataSection + name + "/mode" ).c_str() );
    readComV( dataFile, ( dataSection + name + "/component" ).c_str() );
    readDirection( dataFile, ( dataSection + name + "/direction" ).c_str() );
    readBase3D( dataFile, dataSection + name + "/" );
}

void
BCInterfaceData::showMe( std::ostream& output ) const
{
    output << "baseString = " << M_baseString << std::endl;

    output << "Side       = " << M_side << std::endl;
    output << "Line       = " << M_line << std::endl;
    output << "Quantity   = " << M_quantity << std::endl;
    output << "base       = " << M_base1D.second << std::endl;
    output << "Resistance  = ";
    for ( UInt i(0); i < static_cast<UInt>( M_resistance.size() ); ++i )
        output << M_resistance[i] << " ";
    output << "\n";

    output << "Flag       = " << static_cast< Real > ( M_flag ) << std::endl;
    output << "Type       = " << M_type << std::endl;
    output << "Mode       = " << M_mode << std::endl;
    output << "comV:      = ";
    for ( UInt i(0); i < static_cast<UInt>( M_comV.size() ); ++i )
        output << M_comV[i] << " ";
    output << "\n";
    output << "direction  = " << M_direction << std::endl;
    output << "base       = " << M_base3D.second << std::endl;
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

    M_resistance.clear();
    M_resistance.reserve( resistanceSize );

    for ( UInt j( 0 ); j < resistanceSize; ++j )
        M_resistance.push_back( dataFile( resistance, 0, j ) );
}
void
BCInterfaceData::readBase1D( const GetPot& dataFile, const std::string& base )
{
    for ( std::map< std::string, baseList1D_Type >::iterator j = M_mapBase1D.begin(); j
            != M_mapBase1D.end(); ++j )
        if ( isBase( dataFile, ( base + j->first ).c_str() ) )
        {
            M_base1D.first = j->first;
            M_base1D.second = M_mapBase1D[j->first];

            break;
        }
}

void
BCInterfaceData::readComV( const GetPot& dataFile, const char* component )
{
    UInt componentSize = dataFile.vector_variable_size( component );

    M_comV.clear();
    M_comV.reserve( componentSize );

    for ( UInt j( 0 ); j < componentSize; ++j )
        M_comV.push_back( dataFile( component, 0, j ) );
}

void
BCInterfaceData::readBase3D( const GetPot& dataFile, const std::string& base )
{
    for ( std::map< std::string, baseList3D_Type >::iterator j = M_mapBase3D.begin(); j
            != M_mapBase3D.end(); ++j )
        if ( isBase( dataFile, ( base + j->first ).c_str() ) )
        {
            M_base3D.first = j->first;
            M_base3D.second = M_mapBase3D[j->first];

            break;
        }
}

bool
BCInterfaceData::isBase( const GetPot& dataFile, const char* base )
{
    M_baseString = dataFile( base, " " );

    return dataFile.checkVariable( base );
}

} // Namespace LifeV
