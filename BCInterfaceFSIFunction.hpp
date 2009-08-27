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
   \file BCInterfaceFSIFunction.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-07-15
 */

#ifndef __BCInterfaceFSIFunction_H
#define __BCInterfaceFSIFunction_H 1





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/life.hpp>
#include <life/lifefem/bcFunction.hpp>

#include <life/lifesolver/FSIOperator.hpp>

#include <string>

#include <lifemc/lifefem/BCInterfaceData.hpp>
#include <lifemc/lifefem/BCInterfaceOperatorFunction.hpp>





// ===================================================
//! Namespaces & Enums
// ===================================================
namespace LifeV {





/*!
 * \class BCInterfaceFSIFunction
 * \brief LifeV bcFunction wrapper for BCInterface (FSI problems).
 *
 *  @author Cristiano Malossi
 *  @see
 *
 *  This class is a specialization of BCInterfaceOperatorFunction class for FSI problems.
 *
 *	<b>AVAILABLE OPERATORS</b>
 *
 *	Available operators are:
 *
 *	f_area
 *	f_flux
 *	f_pressure
 *	s_density
 *	s_poisson
 *	s_thickness
 *	s_young
 *
 */
template <class Operator>
class BCInterfaceFSIFunction : 	public virtual BCInterfaceOperatorFunction<Operator>
//     :
//     public LifeV::Application
{
public:

	typedef BCInterfaceOperatorFunction<Operator>		super;

	/** @name Constructors & Destructor
     */
    //@{

    //! Constructor
	BCInterfaceFSIFunction();

    //! Constructor
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	BCInterfaceFSIFunction( const BCInterfaceData<Operator>& data );

	//! Copy constructor
	/*!
	 * \param function			- BCInterfaceFSIFunction
	 */
	BCInterfaceFSIFunction( const BCInterfaceFSIFunction& function );

    //! Destructor
	virtual ~BCInterfaceFSIFunction() {}

    //@}

private:

	/** @name Private functions
	*/
	//@{

	inline void createAccessList( void );

	inline void addOperatorVariables( const Real& t );

	//@}
};

//! Factory create function
template <>
inline BCInterfaceFunction<FSIOperator>* createOperatorFunction<FSIOperator>()
{
	return new BCInterfaceFSIFunction<FSIOperator>();
}



// ===================================================
//! Constructors
// ===================================================
template <class Operator>
BCInterfaceFSIFunction<Operator>::BCInterfaceFSIFunction( ) :
	BCInterfaceFunction<Operator>			( ),
	BCInterfaceOperatorFunction<Operator>	( )
{

#ifdef DEBUG
	Debug( 5025 ) << "BCInterfaceFSIFunction::BCInterfaceFSIFunction( void )" << "\n";
#endif

}



template <class Operator>
BCInterfaceFSIFunction<Operator>::BCInterfaceFSIFunction( const BCInterfaceData<Operator>& data ) :
	BCInterfaceFunction<Operator>			( ),
	BCInterfaceOperatorFunction<Operator>	( )
{

#ifdef DEBUG
	Debug( 5025 ) << "BCInterfaceFSIFunction::BCInterfaceFSIFunction( data )" << "\n";
#endif

	this->setData( data );
}



template <class Operator>
BCInterfaceFSIFunction<Operator>::BCInterfaceFSIFunction( const BCInterfaceFSIFunction& function ) :
	BCInterfaceFunction<Operator>			( function ),
	BCInterfaceOperatorFunction<Operator>	( function )
{
}





// ===================================================
//! Private functions
// ===================================================
template <class Operator>
inline void
BCInterfaceFSIFunction<Operator>::createAccessList( void )
{
#ifdef DEBUG
	Debug( 5025 ) << "BCInterfaceFSIFunction::createAccessList" << "\n";
#endif
	//Create mapList
	super::M_mapList["f_area"]		= super::f_area;
	super::M_mapList["f_flux"]		= super::f_flux;
	super::M_mapList["f_pressure"]	= super::f_pressure;
	super::M_mapList["s_density"]	= super::s_density;
	super::M_mapList["s_poisson"]	= super::s_poisson;
	super::M_mapList["s_thickness"]	= super::s_thickness;
	super::M_mapList["s_young"]		= super::s_young;

	//Create list
	super::createAccessList();
}



template <class Operator>
inline void
BCInterfaceFSIFunction<Operator>::addOperatorVariables( const Real& t )
{

#ifdef DEBUG
	Debug( 5025 ) << "BCInterfaceFSIFunction::addOperatorVariables  " << "\n";
#endif

	//Check if the variables have been already updated
	if ( t == super::M_oldTime )
		return;

	super::M_oldTime = t;

	// Create/Update variables for FSI problem
	for ( typename std::set<typename super::operatorList>::iterator j = super::M_list.begin() ; j != super::M_list.end() ; ++j )
		switch ( *j )
		{
			// f_ -> FLUID
			case super::f_area :

#ifdef DEBUG
				Debug( 5025 ) << "                                                   f_area(" << static_cast<Real> (super::M_flag) << "): " << super::M_operator->fluid().area( super::M_flag )  << "\n";
#endif
				setVariable( "f_area", super::M_operator->fluid().area( super::M_flag ) );

				break;

			case super::f_flux :

#ifdef DEBUG
				Debug( 5025 ) << "                                                   f_flux(" << static_cast<Real> (super::M_flag) << "): " << super::M_operator->fluid().flux( super::M_flag )  << "\n";
#endif

				setVariable( "f_flux", super::M_operator->fluid().flux( super::M_flag ) );

				break;

			case super::f_pressure :

#ifdef DEBUG
				Debug( 5025 ) << "                                               f_pressure(" << static_cast<Real> (super::M_flag) << "): " << super::M_operator->fluid().pressure( super::M_flag )  << "\n";
#endif

				setVariable( "f_pressure", super::M_operator->fluid().pressure( super::M_flag ) );

				break;



			// s_ -> SOLID
			case super::s_density :

#ifdef DEBUG
				Debug( 5025 ) << "                                                   s_density: " << super::M_operator->solid().rho()  << "\n";
#endif

				setVariable( "s_density", super::M_operator->solid().rho() );

				break;

			case super::s_poisson :

#ifdef DEBUG
				Debug( 5025 ) << "                                                   s_poisson: " << super::M_operator->solid().poisson()  << "\n";
#endif

				setVariable( "s_poisson", super::M_operator->solid().poisson() );

				break;

			case super::s_thickness :

#ifdef DEBUG
				Debug( 5025 ) << "                                                 s_thickness: " << super::M_operator->solid().thickness()  << "\n";
#endif

				setVariable( "s_thickness", super::M_operator->solid().thickness() );

				break;

			case super::s_young :

#ifdef DEBUG
				Debug( 5025 ) << "                                                     s_young: " << super::M_operator->solid().young()  << "\n";
#endif

				setVariable( "s_young", super::M_operator->solid().young() );

				break;
		}
}

} // Namespace LifeV

#endif /* __BCInterfaceFSIFunction_H */
