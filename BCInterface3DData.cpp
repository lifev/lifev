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
 *  @brief File containing the BCInterface_Data class
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
BCInterface_Data::BCInterface_Data() :
        M_name               (),
        M_flag               (),
        M_type               (),
        M_mode               (),
        M_comV               (),
        M_direction          (),
        M_baseString         (),
        M_base               (),
        M_mapType            (),
        M_mapMode            (),
        M_mapBase            ()
{
    //Set mapType
    M_mapType["Essential"]  = Essential;
    M_mapType["Natural"]    = Natural;
    M_mapType["Mixte"]      = Mixte;
    M_mapType["Flux"]       = Flux;
    M_mapType["Resistance"] = Resistance;

    //Set mapMode
    M_mapMode["Scalar"]      = Scalar;
    M_mapMode["Full"]        = Full;
    M_mapMode["Component"]   = Component;
    M_mapMode["Normal"]      = Normal;
    M_mapMode["Tangential"]  = Tangential;
    M_mapMode["Directional"] = Directional;

    //Set mapBase
    M_mapBase["function"]         = BCInterface_function;
    M_mapBase["functionFile"]     = BCInterface_functionFile;
    M_mapBase["OPERfunction"]     = BCInterface_OPERfunction;
    M_mapBase["OPERfunctionFile"] = BCInterface_OPERfunctionFile;
    M_mapBase["FSI"]              = BCInterface_OPERFSI;
}

BCInterface_Data::BCInterface_Data( const BCInterface_Data& data ) :
        M_name              ( data.M_name ),
        M_flag              ( data.M_flag ),
        M_type              ( data.M_type ),
        M_mode              ( data.M_mode ),
        M_comV              ( data.M_comV ),
        M_direction         ( data.M_direction ),
        M_baseString        ( data.M_baseString ),
        M_base              ( data.M_base ),
        M_mapType           ( data.M_mapType ),
        M_mapMode           ( data.M_mapMode ),
        M_mapBase           ( data.M_mapBase )
{
}

// ===================================================
// Operators
// ===================================================
BCInterface_Data&
BCInterface_Data::operator=( const BCInterface_Data& data )
{
    if ( this != &data )
    {
        M_name              = data.M_name;
        M_flag              = data.M_flag;
        M_type              = data.M_type;
        M_mode              = data.M_mode;
        M_comV              = data.M_comV;
        M_direction         = data.M_direction;
        M_baseString        = data.M_baseString;
        M_base              = data.M_base;
        M_mapType           = data.M_mapType;
        M_mapMode           = data.M_mapMode;
        M_mapBase           = data.M_mapBase;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
BCInterface_Data::readBC( const std::string& fileName,
                          const std::string& dataSection,
                          const BCName& name )
{
    GetPot dataFile( fileName );

    M_name = name;

    readFlag( dataFile, ( dataSection + name + "/flag" ).c_str() );
    readType( dataFile, ( dataSection + name + "/type" ).c_str() );
    readMode( dataFile, ( dataSection + name + "/mode" ).c_str() );
    readComV( dataFile, ( dataSection + name + "/component" ).c_str() );
    readDirection( dataFile, ( dataSection + name + "/direction" ).c_str() );
    readBase( dataFile, dataSection + name + "/" );
}

void
BCInterface_Data::showMe( std::ostream& output ) const
{
    output << "Flag       = " << static_cast< Real > ( M_flag ) << std::endl;
    output << "Type       = " << M_type << std::endl;
    output << "Mode       = " << M_mode << std::endl;
    output << "comV:      = ";
    for ( UInt i(0); i < static_cast<UInt>( M_comV.size() ); ++i )
        output << M_comV[i] << " ";
    output << "\n";
    output << "direction  = " << M_direction << std::endl;
    output << "base       = " << M_base.second << std::endl;
    output << "baseString = " << M_baseString << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void BCInterface_Data::setBaseString( const std::string& baseString )
{
    M_baseString = baseString;
    boost::replace_all( M_baseString, " ", "" );
}

// ===================================================
// Private Methods
// ===================================================
void
BCInterface_Data::readComV( const GetPot& dataFile, const char* component )
{
    UInt componentSize = dataFile.vector_variable_size( component );

    M_comV.clear();
    M_comV.reserve( componentSize );

    for ( UInt j( 0 ); j < componentSize; ++j )
        M_comV.push_back( dataFile( component, 0, j ) );
}

void
BCInterface_Data::readBase( const GetPot& dataFile, const std::string& base )
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
BCInterface_Data::isBase( const GetPot& dataFile, const char* base )
{
    M_baseString = dataFile( base, " " );

    return dataFile.checkVariable( base );
}

} // Namespace LifeV
