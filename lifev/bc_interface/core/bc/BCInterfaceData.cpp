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

#include <lifev/bc_interface/core/bc/BCInterfaceData.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterfaceData::BCInterfaceData() :
    M_base                  (),
    M_baseString            (),
    M_mapBase               (),
    M_boundaryID            (),
    M_parameters            ()
{
    //Set mapBase
    M_mapBase["function"]           = BCIFunctionParser;
    M_mapBase["functionFile"]       = BCIFunctionParserFile;
    M_mapBase["functionSolver"]     = BCIFunctionParserSolver;
    M_mapBase["functionFileSolver"] = BCIFunctionParserFileSolver;
    M_mapBase["functionUD"]         = BCIFunctionUserDefined;
    M_mapBase["functionSD"]         = BCIFunctionSolverDefined;
    M_mapBase["dataInterpolator"]   = BCI3DDataInterpolator;
}

BCInterfaceData::BCInterfaceData ( const BCInterfaceData& data ) :
    M_base                  ( data.M_base ),
    M_baseString            ( data.M_baseString ),
    M_mapBase               ( data.M_mapBase ),
    M_boundaryID            ( data.M_boundaryID ),
    M_parameters            ( data.M_parameters )
{
}

// ===================================================
// Operators
// ===================================================
BCInterfaceData&
BCInterfaceData::operator= ( const BCInterfaceData& data )
{
    if ( this != &data )
    {
        M_base       = data.M_base;
        M_baseString = data.M_baseString;
        M_mapBase    = data.M_mapBase;
        M_boundaryID = data.M_boundaryID;
        M_parameters = data.M_parameters;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
BCInterfaceData::readBC ( const std::string& fileName, const std::string& dataSection, const std::string& name )
{
    GetPot dataFile ( fileName );

    // Read base
    readBase ( dataFile, dataSection + name + "/", M_base, M_baseString );

    // Read parameters
    readParameters ( dataFile, ( dataSection + name + "/parameters" ).c_str() );
}

void
BCInterfaceData::showMe ( std::ostream& output ) const
{
    output << "baseString        = " << M_baseString << std::endl;
    output << "base              = " << M_base.second << std::endl;

    output << "boundary ID       = " << M_boundaryID << std::endl;

    output << "Parameters        = ";
    for ( UInt i (0); i < static_cast<UInt> ( M_parameters.size() ); ++i )
    {
        output << M_parameters[i] << " ";
    }
    output << "\n";
}

// ===================================================
// Set Methods
// ===================================================
void BCInterfaceData::setBaseString ( const std::string& baseString )
{
    M_baseString = baseString;
	std::string search = " ";
	std::string replace = "";

    for( size_t pos = 0; ; pos += replace.length() ) {
        pos = M_baseString.find( search, pos );
        if( pos == std::string::npos ) break;
        M_baseString.erase( pos, search.length() );
        M_baseString.insert( pos, replace );
    }
}

// ===================================================
// Private Methods
// ===================================================
void
BCInterfaceData::readBase ( const GetPot& dataFile, const std::string& path, std::pair< std::string, baseList_Type >& base, std::string& baseString )
{
    for ( std::map< std::string, baseList_Type >::iterator j = M_mapBase.begin(); j != M_mapBase.end(); ++j )
        if ( isBase ( dataFile, ( path + j->first ).c_str(), baseString ) )
        {
            base.first = j->first;
            base.second = M_mapBase[j->first];

            break;
        }
}

bool
BCInterfaceData::isBase ( const GetPot& dataFile, const char* base, std::string& baseString )
{
    baseString = dataFile ( base, " " );

    return dataFile.checkVariable ( base );
}

void
BCInterfaceData::readParameters ( const GetPot& dataFile, const char* parameters )
{
    UInt parametersSize = dataFile.vector_variable_size ( parameters );

    M_parameters.resize ( parametersSize );
    for ( UInt j ( 0 ); j < parametersSize; ++j )
    {
        M_parameters[j] = dataFile ( parameters, 0, j );
    }
}

} // Namespace LifeV
