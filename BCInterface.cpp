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

// Initialize static variables
std::vector< boost::shared_ptr<BCInterfaceFunction> >			BCInterface::M_vectorFunction;
std::map<std::string,size_type>									BCInterface::M_mapFunction;

std::vector< boost::shared_ptr<BCInterfaceFunctionFile> >		BCInterface::M_vectorFunctionFile;
std::map<std::string,size_type>									BCInterface::M_mapFunctionFile;

std::vector< boost::shared_ptr<BCInterfaceFSIFunction> >		BCInterface::M_vectorFSIFunction;
std::map<std::string,size_type>									BCInterface::M_mapFSIFunction;

std::vector< boost::shared_ptr<BCInterfaceFSIFunctionFile> >	BCInterface::M_vectorFSIFunctionFile;
std::map<std::string,size_type>									BCInterface::M_mapFSIFunctionFile;

// ===================================================
//! Constructor
// ===================================================
BCInterface::BCInterface( const GetPot& dataFile, const std::string& dataSection ) :
	M_dataFile					( dataFile ),
	M_dataSection				( dataSection + "/boundary_conditions/" ),
	M_list						( ),
	M_listSize					( 0 ),
	M_autoSetParameters			( true ),
	M_bcNumber					( 0 ),
	M_hint						( BCHandler::HINT_BC_ONLY_ESSENTIAL ),
	M_handler					( ),
	M_mapType					( ),
	M_mapMode					( ),
	M_mapBase					( ),
	M_data						( ),
	M_base						( ),
	M_FSIOperator				( ),
	M_vectorFSI					( )

{

#ifdef DEBUG
	    Debug( 5020 ) << "BCInterface::BCInterface------------------------------" << "\n";
#endif

	//Set mapType
	M_mapType["Essential"] 			= Essential;
	M_mapType["Natural"] 			= Natural;
	M_mapType["Mixte"] 				= Mixte;
	//M_mapType["Flux"] 			= Flux;

	//Set mapMode
	M_mapMode["Scalar"] 			= Scalar;
	M_mapMode["Full"] 				= Full;
	M_mapMode["Component"] 			= Component;
	M_mapMode["Normal"] 			= Normal;
	M_mapMode["Tangential"] 		= Tangential;

	//Set mapBase
	M_mapBase["function"] 			= function;
	M_mapBase["functionFile"] 		= functionFile;
	M_mapBase["FSI"]				= FSI;
	M_mapBase["FSIfunction"]		= FSIfunction;
	M_mapBase["FSIfunctionFile"]	= FSIfunctionFile;

	//Clear non-static vectors
	M_vectorFSI.clear();

	//Set other parameters
	setList( (M_dataSection + "list").c_str() );
}


BCInterface::BCInterface( const BCInterface& interface ) :
	M_dataFile				( interface.M_dataFile ),
	M_dataSection			( interface.M_dataSection ),
	M_list					( interface.M_list ),
	M_listSize				( interface.M_listSize ),
	M_autoSetParameters		( interface.M_autoSetParameters ),
	M_bcNumber				( interface.M_bcNumber ),
	M_hint					( interface.M_hint ),
	M_handler				( interface.M_handler ),
	M_mapType				( interface.M_mapType ),
	M_mapMode				( interface.M_mapMode ),
	M_mapBase				( interface.M_mapBase ),
	M_data					( interface.M_data ),
	M_base					( interface.M_base ),
	M_FSIOperator			( interface.M_FSIOperator ),
	M_vectorFSI				( interface.M_vectorFSI )
{
}



BCInterface&
BCInterface::operator=( const BCInterface& interface )
{
    if ( this != &interface )
    {
    	M_dataFile				= interface.M_dataFile;
    	M_dataSection			= interface.M_dataSection;
    	M_list					= interface.M_list;
    	M_listSize				= interface.M_listSize;
    	M_autoSetParameters		= interface.M_autoSetParameters;
    	M_bcNumber				= interface.M_bcNumber;
    	M_hint					= interface.M_hint;
    	M_handler				= interface.M_handler;
    	M_mapType				= interface.M_mapType;
    	M_mapMode				= interface.M_mapMode;
    	M_mapBase				= interface.M_mapBase;
    	M_data					= interface.M_data;
    	M_base					= interface.M_base;
    	M_FSIOperator			= interface.M_FSIOperator;
    	M_vectorFSI				= interface.M_vectorFSI;
    }

	return *this;
}





// ===================================================
//! Members functions
// ===================================================
void
BCInterface::setHandlerParameters( const ID& bcNumber, const BCHandler::BCHints& hint )
{
	M_bcNumber 	= bcNumber;
	M_hint 		= hint;

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface::setHandlerParameters          M_bcNumber: " << M_bcNumber << "\n";
    Debug( 5020 ) << "                                               M_hint: " << M_hint << "\n";
#endif

	M_autoSetParameters = false;
}



void
BCInterface::buildHandler( void )
{

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::buildHandler         M_autoSetParameters: " << M_autoSetParameters << "\n";
#endif

	if ( M_autoSetParameters )
		autosetHandlerParameters();

	M_handler.reset( new BCHandler ( M_bcNumber, M_hint ) );

	for ( UInt i(0) ; i < M_listSize ; ++i )
	{
		M_data.set_name( M_list[i] );

		readFlag( (M_dataSection + M_data.get_name() + "/flag").c_str() );
		readType( (M_dataSection + M_data.get_name() + "/type").c_str() );
		readMode( (M_dataSection + M_data.get_name() + "/mode").c_str() );
		readComV( (M_dataSection + M_data.get_name() + "/component").c_str() );
		readBase(  M_dataSection + M_data.get_name() + "/" );

		switch ( M_base )
		{
			case function :

				if ( newBase( M_mapFunction, M_vectorFunction ) )
					addBase( M_vectorFunction );

				addBCManager( M_vectorFunction[M_mapFunction[M_data.get_baseString()]]->getBase() );

				break;

			case functionFile :

				if ( newBase( M_mapFunctionFile, M_vectorFunctionFile ) )
					addBase( M_vectorFunctionFile );

				addBCManager( M_vectorFunctionFile[M_mapFunctionFile[M_data.get_baseString()]]->getBase() );

				break;

			case FSI :

				addBase( M_vectorFSI, M_FSIOperator );
				addBCManager( M_vectorFSI.back()->getBase() );

				break;

			case FSIfunction :

				if ( newBase( M_mapFSIFunction, M_vectorFSIFunction ) )
					addBase( M_vectorFSIFunction, M_FSIOperator );

				addBCManager( M_vectorFSIFunction.back()->getBase() );

				break;

			case FSIfunctionFile :

				if ( newBase( M_mapFSIFunctionFile, M_vectorFSIFunctionFile ) )
					addBase( M_vectorFSIFunctionFile, M_FSIOperator );

				addBCManager( M_vectorFSIFunctionFile.back()->getBase() );

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
		if ( M_data.get_type() != Essential )
			M_hint = BCHandler::HINT_BC_NONE;

		M_bcNumber += M_dataFile.vector_variable_size((M_dataSection + M_list[i] + "/flag").c_str());
	}

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface::autosetHandlerParameters      M_bcNumber: " << M_bcNumber << "\n";
    Debug( 5020 ) << "                                               M_hint: " << M_hint << "\n\n";
#endif

}



inline void
BCInterface::readFlag( const char* flag )
{
	M_data.set_flag( M_dataFile(flag, 0) );

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface::readFlag                            flag: " << static_cast<Real>(M_data.get_flag()) << "\n";
#endif
}



inline void
BCInterface::readType( const char* type )
{
	M_data.set_type( M_mapType[M_dataFile(type, "Essential")] );

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::readType                            type: " << M_data.get_type() << " (" << M_dataFile(type, "Essential") << ")\n";
#endif
}



inline void
BCInterface::readMode( const char* mode )
{
	M_data.set_mode( M_mapMode[M_dataFile(mode, "Full")] );

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::readMode                            mode: " << M_data.get_mode() << " (" << M_dataFile(mode, "Full") << ")\n";
#endif
}


inline void
BCInterface::readComV( const char* component )
{
    UInt componentSize = M_dataFile.vector_variable_size(component);

    M_data.reset_comV( componentSize );

    for (UInt j(0) ; j < componentSize ; ++j)
    	M_data.set_comV( M_dataFile(component, 0, j) );

#ifdef DEBUG
    std::stringstream output;
    output << "BCInterface::readComV                            comV: ";
    for (UInt i(0) ; i < static_cast<UInt>(M_data.get_comV().size()) ; ++i )
    	output << M_data.get_comV()[i] << " ";
    Debug( 5020 ) << output.str() << "\n";
#endif

}



inline void
BCInterface::readBase( const std::string& base )
{
	for ( std::map<std::string, BCBaseList>::iterator j = M_mapBase.begin() ; j != M_mapBase.end() ; ++j )
		if ( isBase( (base + j->first).c_str() ) )
		{
			M_base = M_mapBase[j->first];

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::readBase                            base: " << M_base << " (" << j->first << ")\n";
			Debug( 5020 ) << "                                           baseString: " << M_data.get_baseString() << "\n";
#endif

			break;
		}
}



inline bool
BCInterface::isBase( const char* base )
{
	M_data.set_baseString( M_dataFile( base, " " ) );

	return M_dataFile.checkVariable( base );
}

} // Namespace LifeV
