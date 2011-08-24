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
        M_mapBase               ()
{
    //Set mapBase
    M_mapBase["function"]           = BCIFunctionParser;
    M_mapBase["functionFile"]       = BCIFunctionParserFile;
    M_mapBase["functionSolver"]     = BCIFunctionParserSolver;
    M_mapBase["functionFileSolver"] = BCIFunctionParserFileSolver;
    M_mapBase["functionUD"]         = BCIFunctionUserDefined;
    M_mapBase["functionDefault"]    = BCI1DFunctionDefault;
    M_mapBase["dataInterpolator"]   = BCI3DDataInterpolator;
}

BCInterfaceData::BCInterfaceData( const BCInterfaceData& data ) :
        M_base                  ( data.M_base ),
        M_baseString            ( data.M_baseString ),
        M_mapBase               ( data.M_mapBase )
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
        M_mapBase               = data.M_mapBase;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
BCInterfaceData::readBC( const std::string& fileName, const std::string& dataSection, const std::string& name )
{
    GetPot dataFile( fileName );

    // Read base
    readBase( dataFile, dataSection + name + "/", M_base, M_baseString );
}

void
BCInterfaceData::showMe( std::ostream& output ) const
{
    output << "baseString        = " << M_baseString << std::endl;
    output << "base              = " << M_base.second << std::endl;
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
