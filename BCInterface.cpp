/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-04-01

  Copyright (C) 2009 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file BCInterface.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-04-01
 */

#include <lifemc/lifefem/BCInterface.hpp>

namespace LifeV {

// ===================================================
//! Constructor
// ===================================================
BCInterface::BCInterface( GetPot const& dataFile, const std::string dataSection ) :
	M_dataFile					( dataFile ),
	M_dataSection				( dataSection + "/boundary_conditions/" ),
	M_list						( ),
	M_listSize					( 0 ),
	M_autoSetParameters			( true ),
	M_bcNumber					( 0 ),
	M_hint						( BCHandler::HINT_BC_ONLY_ESSENTIAL ),
	M_handler					( ),
	M_FSIOperator				( ),
	M_name						( "undefined" ),
	M_flag						( ),
	M_type						( ),
	M_mode						( ),
	M_comN						( ),
	M_comV						( ),
	M_base						( ),
	M_baseString				( "undefined" )
{
	//Set mapType
	M_mapType["Essential"] 	= Essential;
	M_mapType["Natural"] 	= Natural;
	M_mapType["Mixte"] 		= Mixte;
	//M_mapType["Flux"] 		= Flux;

	//Set mapMode
	M_mapMode["Scalar"] 	= Scalar;
	M_mapMode["Full"] 		= Full;
	M_mapMode["Component"] 	= Component;
	M_mapMode["Normal"] 	= Normal;
	M_mapMode["Tangential"] = Tangential;

	//Set mapBase
	M_mapBase["function"] 	= function;
	M_mapBase["fsi"] 		= fsi;

	//Operators
	M_functionVector.clear();
	M_FSIOperatorVector.clear();

	//Set other parameters
	setList( (M_dataSection + "list").c_str() );
}





// ===================================================
//! Members functions
// ===================================================
void
BCInterface::setHandlerParameters( const ID bcNumber, const BCHandler::BCHints hint )
{
	M_bcNumber 	= bcNumber;
	M_hint 		= hint;

    Debug( 5020 ) << "BCInterface::setHandlerParameters          M_bcNumber: " << M_bcNumber << "\n";
    Debug( 5020 ) << "                                               M_hint: " << M_hint << "\n";

	M_autoSetParameters = false;
}



void
BCInterface::setFSIOperator( const boost::shared_ptr<FSIOperator>& oper )
{
	M_FSIOperator = oper;
}



void
BCInterface::buildHandler( void )
{
	Debug( 5020 ) << "BCInterface::buildHandler         M_autoSetParameters: " << M_autoSetParameters << "\n";

	if ( M_autoSetParameters )
		autosetHandlerParameters();

	M_handler.reset( new BCHandler ( M_bcNumber, M_hint ) );

	for ( UInt i(0) ; i < M_listSize ; ++i )
	{
		M_name = M_list[i];

		readFlag( (M_dataSection + M_name + "/flag").c_str() );
		readType( (M_dataSection + M_name + "/type").c_str() );
		readMode( (M_dataSection + M_name + "/mode").c_str() );

		readBase(  M_dataSection + M_name + "/" );

		switch ( M_base )
		{
			case function :

				addBase( M_functionVector );
				addBCManager( M_functionVector.back()->getBase() );

				break;

			case fsi :

				addBase( M_FSIOperatorVector, M_FSIOperator );
				addBCManager( M_FSIOperatorVector.back()->getBase() );

				break;
		}
	}
}





// ===================================================
//! Private functions
// ===================================================
inline void
BCInterface::setList( const char* conditions )
{
    M_listSize = M_dataFile.vector_variable_size( conditions );

    M_list.reserve( M_listSize );
    for ( UInt i(0) ; i < M_listSize ; ++i )
    	M_list.push_back(M_dataFile(conditions, " ", i));
}



inline void
BCInterface::autosetHandlerParameters( void )
{
	for ( UInt i(0) ; i < M_listSize ; ++i )
	{
		readType( (M_dataSection + M_list[i] + "/type").c_str() );
		if ( M_type != Essential )
			M_hint = BCHandler::HINT_BC_NONE;

		M_bcNumber += M_dataFile.vector_variable_size((M_dataSection + M_list[i] + "/flag").c_str());
	}

    Debug( 5020 ) << "BCInterface::autosetHandlerParameters      M_bcNumber: " << M_bcNumber << "\n";
    Debug( 5020 ) << "                                               M_hint: " << M_hint << "\n\n";
}



inline void
BCInterface::readFlag( const char* flag )
{
    UInt flagSize = M_dataFile.vector_variable_size(flag);

	M_flag.clear();
	M_flag.reserve( flagSize );

    for ( UInt j(0) ; j < flagSize ; ++j )
    	M_flag.push_back( M_dataFile(flag, 0, j) );

    Debug( 5020 ) << "BCInterface::readFlag                   M_flag.size(): " << Real(M_flag.size()) << "\n";
}



inline void
BCInterface::readType( const char* type )
{
	M_type = M_mapType[M_dataFile(type, "Essential")];

	Debug( 5020 ) << "BCInterface::readType                          M_type: " << M_type << " " << M_dataFile(type, "Essential") << "\n";
}



inline void
BCInterface::readMode( const char* mode )
{
	M_mode = M_mapMode[M_dataFile(mode, "Full")];

	Debug( 5020 ) << "BCInterface::readMode                          M_mode: " << M_mode << " " << M_dataFile(mode, "Full") << "\n";
}



inline void
BCInterface::readComponentNumber( const char* component )
{
	M_comN = M_dataFile( component, 0 );

	Debug( 5020 ) << "BCInterface::readComponentNumber               M_comN: " << M_comN << "\n";
}



inline void
BCInterface::readComponentVector( const char* component )
{
    UInt componentSize = M_dataFile.vector_variable_size(component);

	M_comV.clear();
    M_comV.reserve( componentSize );

    for (UInt j(0) ; j < componentSize ; ++j)
    	M_comV.push_back( M_dataFile(component, 0, j) );

    Debug( 5020 ) << "BCInterface::readComponentVector        M_comV.size(): " << Real(M_comV.size()) << "\n";
}



inline void
BCInterface::readBase( const std::string base )
{
	for ( std::map<std::string, BCBaseList>::iterator j = M_mapBase.begin() ; j != M_mapBase.end() ; ++j )
		if ( isBase( (base + j->first).c_str() ) )
		{
			M_base = M_mapBase[j->first];

			Debug( 5020 ) << "BCInterface::readBase                          M_base: " << M_base << " " << j->first << "\n";
			Debug( 5020 ) << "                                         M_baseString: " << M_baseString << "\n";

			break;
		}
}



inline bool
BCInterface::isBase( const char* base )
{
	M_baseString = M_dataFile( base, " " );

	return M_dataFile.checkVariable( base );
}

} // Namespace LifeV
