/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-07-15

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
   \file BCInterfaceFSIFunction.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-07-15
 */

#include <lifemc/lifefem/BCInterfaceFSIFunction.hpp>

namespace LifeV {

// ===================================================
//! Constructor & Destructor
// ===================================================
BCInterfaceFSIFunction::BCInterfaceFSIFunction( const boost::shared_ptr<FSIOperator>& oper ) :
	BCInterfaceFunction		( ),
	M_FSIOperator			( oper ),
	M_flag					( ),
	M_mapFSIList			( ),
	M_FSIList				( ),
	M_oldTime				( -1.0 ) // Negative time
{

#ifdef DEBUG
	Debug( 5024 ) << "BCInterfaceFSIFunction::BCInterfaceFSIFunction( void )" << "\n";
#endif

}



BCInterfaceFSIFunction::BCInterfaceFSIFunction( const BCInterfaceData& data,
												const boost::shared_ptr<FSIOperator>& oper ) :
	BCInterfaceFunction		( ),
	M_FSIOperator			( oper ),
	M_flag					( ),
	M_mapFSIList			( ),
	M_FSIList				( ),
	M_oldTime				( -1.0 ) // Negative time
{

#ifdef DEBUG
	Debug( 5024 ) << "BCInterfaceFSIFunction::BCInterfaceFSIFunction" << "\n";
#endif

	this->setData( data );
}



BCInterfaceFSIFunction::BCInterfaceFSIFunction( const BCInterfaceFSIFunction& function ) :
	BCInterfaceFunction		( function ),
	M_FSIOperator			( function.M_FSIOperator ),
	M_flag					( function.M_flag ),
	M_mapFSIList			( function.M_mapFSIList ),
	M_FSIList				( function.M_FSIList ),
	M_oldTime				( function.M_oldTime )
{
}





// ===================================================
//! Methods
// ===================================================
BCInterfaceFSIFunction&
BCInterfaceFSIFunction::operator=( const BCInterfaceFSIFunction& function )
{
    if ( this != &function )
    {
    	BCInterfaceFunction::operator=( function );
    	M_FSIOperator	= function.M_FSIOperator;
    	M_flag			= function.M_flag;
    	M_mapFSIList	= function.M_mapFSIList;
    	M_FSIList		= function.M_FSIList;
    	M_oldTime		= function.M_oldTime;
    }

	return *this;
}



void
BCInterfaceFSIFunction::setData( const BCInterfaceData& data )
{

#ifdef DEBUG
	Debug( 5024 ) << "BCInterfaceFSIFunction::setData" << "\n";
#endif

	M_flag = data.get_flag();

	BCInterfaceFunction::setData( data );

	createFSIAccessList();
}



bool
BCInterfaceFSIFunction::compare( const BCInterfaceData& data )
{
	return M_baseString.compare( data.get_baseString() ) == 0 && M_comV == data.get_comV() && M_flag == data.get_flag();
}





// ===================================================
//! Private functions
// ===================================================
inline void
BCInterfaceFSIFunction::createFSIAccessList( void )
{
	//Create mapFSIList
	M_mapFSIList["f_area"]		= f_area;
	M_mapFSIList["f_flux"]		= f_flux;
	M_mapFSIList["f_pressure"]	= f_pressure;
	M_mapFSIList["s_density"]	= s_density;
	M_mapFSIList["s_poisson"]	= s_poisson;
	M_mapFSIList["s_thickness"]	= s_thickness;
	M_mapFSIList["s_young"]		= s_young;

	//Create FSIList
	M_FSIList.clear();
	for ( std::map<std::string, FSIList>::iterator j = M_mapFSIList.begin() ; j != M_mapFSIList.end() ; ++j )
		if ( boost::find_first( M_baseString, j->first ) )
			M_FSIList.insert( j->second );
}



inline void
BCInterfaceFSIFunction::addFSIVariables( const Real& t )
{

#ifdef DEBUG
	Debug( 5024 ) << "BCInterfaceFSIFunction::addFSIVariables  " << "\n";
#endif

	//Check if FSIVariables have to been updated
	if ( t == M_oldTime )
		return;
	M_oldTime = t;


	// Create/Update FSIVariables
	for ( std::set<FSIList>::iterator j = M_FSIList.begin() ; j != M_FSIList.end() ; ++j )
		switch ( *j )
		{
			// f_ -> FLUID
			case f_area :

#ifdef DEBUG
				Debug( 5024 ) << "                                                   f_area(" << static_cast<Real> (M_flag) << "): " << M_FSIOperator->fluid().area( M_flag )  << "\n";
#endif
				M_parser->setVariable( "f_area", M_FSIOperator->fluid().area( M_flag ) );

				break;

			case f_flux :

#ifdef DEBUG
				Debug( 5024 ) << "                                                   f_flux(" << static_cast<Real> (M_flag) << "): " << M_FSIOperator->fluid().flux( M_flag )  << "\n";
#endif

				M_parser->setVariable( "f_flux", M_FSIOperator->fluid().flux( M_flag ) );

				break;

			case f_pressure :

#ifdef DEBUG
				Debug( 5024 ) << "                                               f_pressure(" << static_cast<Real> (M_flag) << "): " << M_FSIOperator->fluid().pressure( M_flag )  << "\n";
#endif

				M_parser->setVariable( "f_pressure", M_FSIOperator->fluid().pressure( M_flag ) );

				break;



			// s_ -> SOLID
			case s_density :

#ifdef DEBUG
				Debug( 5024 ) << "                                                   s_density: " << M_FSIOperator->solid().rho()  << "\n";
#endif

				M_parser->setVariable( "s_density", M_FSIOperator->solid().rho() );

				break;

			case s_poisson :

#ifdef DEBUG
				Debug( 5024 ) << "                                                   s_poisson: " << M_FSIOperator->solid().poisson()  << "\n";
#endif

				M_parser->setVariable( "s_poisson", M_FSIOperator->solid().poisson() );

				break;

			case s_thickness :

#ifdef DEBUG
				Debug( 5024 ) << "                                                 s_thickness: " << M_FSIOperator->solid().thickness()  << "\n";
#endif

				M_parser->setVariable( "s_thickness", M_FSIOperator->solid().thickness() );

				break;

			case s_young :

#ifdef DEBUG
				Debug( 5024 ) << "                                                     s_young: " << M_FSIOperator->solid().young()  << "\n";
#endif

				M_parser->setVariable( "s_young", M_FSIOperator->solid().young() );

				break;
		}
}

} // Namespace LifeV
