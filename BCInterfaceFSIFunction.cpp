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
BCInterfaceFSIFunction::BCInterfaceFSIFunction( ) :
	BCInterfaceFunction							( ),
	BCInterfaceOperatorFunction<FSIOperator>	( )
{

#ifdef DEBUG
	Debug( 5025 ) << "BCInterfaceFSIFunction::BCInterfaceFSIFunction( void )" << "\n";
#endif

}



BCInterfaceFSIFunction::BCInterfaceFSIFunction( const boost::shared_ptr<FSIOperator>& Oper ) :
	BCInterfaceFunction							( ),
	BCInterfaceOperatorFunction<FSIOperator>	( Oper )
{

#ifdef DEBUG
	Debug( 5025 ) << "BCInterfaceFSIFunction::BCInterfaceFSIFunction( Oper )" << "\n";
#endif

}



BCInterfaceFSIFunction::BCInterfaceFSIFunction( const BCInterfaceData& data,
												const boost::shared_ptr<FSIOperator>& Oper ) :
	BCInterfaceFunction							( ),
	BCInterfaceOperatorFunction<FSIOperator>	( Oper )
{

#ifdef DEBUG
	Debug( 5025 ) << "BCInterfaceFSIFunction::BCInterfaceFSIFunction( data, Oper )" << "\n";
#endif

	this->setData( data );
}



BCInterfaceFSIFunction::BCInterfaceFSIFunction( const BCInterfaceFSIFunction& function ) :
	BCInterfaceFunction							( function ),
	BCInterfaceOperatorFunction<FSIOperator>	( function )
{
}





// ===================================================
//! Private functions
// ===================================================
void
BCInterfaceFSIFunction::createAccessList( void )
{
#ifdef DEBUG
	Debug( 5025 ) << "BCInterfaceFSIFunction::createAccessList" << "\n";
#endif
	//Create mapList
	M_mapList["f_area"]			= f_area;
	M_mapList["f_flux"]			= f_flux;
	M_mapList["f_pressure"]		= f_pressure;
	M_mapList["s_density"]		= s_density;
	M_mapList["s_poisson"]		= s_poisson;
	M_mapList["s_thickness"]	= s_thickness;
	M_mapList["s_young"]		= s_young;

	//Create list
	BCInterfaceOperatorFunction<FSIOperator>::createAccessList();
}



void
BCInterfaceFSIFunction::addOperatorVariables( const Real& t )
{

#ifdef DEBUG
	Debug( 5025 ) << "BCInterfaceFSIFunction::addOperatorVariables  " << "\n";
#endif

	BCInterfaceOperatorFunction<FSIOperator>::addOperatorVariables( t );

	// Create/Update variables for FSI problem
	for ( std::set<operatorList>::iterator j = M_list.begin() ; j != M_list.end() ; ++j )
		switch ( *j )
		{
			// f_ -> FLUID
			case f_area :

#ifdef DEBUG
				Debug( 5025 ) << "                                                   f_area(" << static_cast<Real> (M_flag) << "): " << M_operator->fluid().area( M_flag )  << "\n";
#endif
				setVariable( "f_area", M_operator->fluid().area( M_flag ) );

				break;

			case f_flux :

#ifdef DEBUG
				Debug( 5025 ) << "                                                   f_flux(" << static_cast<Real> (M_flag) << "): " << M_operator->fluid().flux( M_flag )  << "\n";
#endif

				setVariable( "f_flux", M_operator->fluid().flux( M_flag ) );

				break;

			case f_pressure :

#ifdef DEBUG
				Debug( 5025 ) << "                                               f_pressure(" << static_cast<Real> (M_flag) << "): " << M_operator->fluid().pressure( M_flag )  << "\n";
#endif

				setVariable( "f_pressure", M_operator->fluid().pressure( M_flag ) );

				break;



			// s_ -> SOLID
			case s_density :

#ifdef DEBUG
				Debug( 5025 ) << "                                                   s_density: " << M_operator->solid().rho()  << "\n";
#endif

				setVariable( "s_density", M_operator->solid().rho() );

				break;

			case s_poisson :

#ifdef DEBUG
				Debug( 5025 ) << "                                                   s_poisson: " << M_operator->solid().poisson()  << "\n";
#endif

				setVariable( "s_poisson", M_operator->solid().poisson() );

				break;

			case s_thickness :

#ifdef DEBUG
				Debug( 5025 ) << "                                                 s_thickness: " << M_operator->solid().thickness()  << "\n";
#endif

				setVariable( "s_thickness", M_operator->solid().thickness() );

				break;

			case s_young :

#ifdef DEBUG
				Debug( 5025 ) << "                                                     s_young: " << M_operator->solid().young()  << "\n";
#endif

				setVariable( "s_young", M_operator->solid().young() );

				break;
		}
}

} // Namespace LifeV
