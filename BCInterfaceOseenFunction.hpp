/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-08-27

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
   \date 2009-08-27
 */

#ifndef __BCInterfaceOseenFunction_H
#define __BCInterfaceOseenFunction_H 1





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/life.hpp>
#include <life/lifefem/bcFunction.hpp>

#include <life/lifesolver/Oseen.hpp>

#include <string>

#include <lifemc/lifefem/BCInterfaceData.hpp>
#include <lifemc/lifefem/BCInterfaceOperatorFunction.hpp>





// ===================================================
//! Namespaces & Enums
// ===================================================
namespace LifeV {





/*!
 * \class BCInterfaceOseenFunction
 * \brief LifeV bcFunction wrapper for BCInterface (Oseen problems).
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
 *
 */
template <class Operator>
class BCInterfaceOseenFunction : 	public virtual BCInterfaceOperatorFunction<Operator>
//     :
//     public LifeV::Application
{
public:

	typedef BCInterfaceOperatorFunction<Operator>		super;

	/** @name Constructors & Destructor
     */
    //@{

    //! Constructor
	BCInterfaceOseenFunction();

    //! Constructor
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	BCInterfaceOseenFunction( const BCInterfaceData<Operator>& data );

	//! Copy constructor
	/*!
	 * \param function			- BCInterfaceOseenFunction
	 */
	BCInterfaceOseenFunction( const BCInterfaceOseenFunction& function );

    //! Destructor
	virtual ~BCInterfaceOseenFunction() {}

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
template <> template <typename Mesh, typename SolverType>
inline BCInterfaceFunction< Oseen<Mesh, SolverType> >* createOperatorFunction< Oseen<Mesh, SolverType> >()
{
	return new BCInterfaceOseenFunction< Oseen<Mesh, SolverType> >();
}



// ===================================================
//! Constructors
// ===================================================
template <class Operator>
BCInterfaceOseenFunction<Operator>::BCInterfaceOseenFunction( ) :
	BCInterfaceFunction<Operator>			( ),
	BCInterfaceOperatorFunction<Operator>	( )
{

#ifdef DEBUG
	Debug( 5027 ) << "BCInterfaceOseenFunction::BCInterfaceOseenFunction( void )" << "\n";
#endif

}



template <class Operator>
BCInterfaceOseenFunction<Operator>::BCInterfaceOseenFunction( const BCInterfaceData<Operator>& data ) :
	BCInterfaceFunction<Operator>			( ),
	BCInterfaceOperatorFunction<Operator>	( )
{

#ifdef DEBUG
	Debug( 5027 ) << "BCInterfaceOseenFunction::BCInterfaceOseenFunction( data )" << "\n";
#endif

	this->setData( data );
}



template <class Operator>
BCInterfaceOseenFunction<Operator>::BCInterfaceOseenFunction( const BCInterfaceOseenFunction& function ) :
	BCInterfaceFunction<Operator>			( function ),
	BCInterfaceOperatorFunction<Operator>	( function )
{
}





// ===================================================
//! Private functions
// ===================================================
template <class Operator>
inline void
BCInterfaceOseenFunction<Operator>::createAccessList( void )
{
#ifdef DEBUG
	Debug( 5027 ) << "BCInterfaceOseenFunction::createAccessList" << "\n";
#endif
	//Create mapList
	super::M_mapList["f_area"]		= super::f_area;
	super::M_mapList["f_flux"]		= super::f_flux;
	super::M_mapList["f_pressure"]	= super::f_pressure;

	//Create list
	super::createAccessList();
}



template <class Operator>
inline void
BCInterfaceOseenFunction<Operator>::addOperatorVariables( const Real& t )
{

#ifdef DEBUG
	Debug( 5027 ) << "BCInterfaceOseenFunction::addOperatorVariables  " << "\n";
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
				Debug( 5027 ) << "                                                   f_area(" << static_cast<Real> (super::M_flag) << "): " << super::M_operator->area( super::M_flag )  << "\n";
#endif
				setVariable( "f_area", super::M_operator->area( super::M_flag ) );

				break;

			case super::f_flux :

#ifdef DEBUG
				Debug( 5027 ) << "                                                   f_flux(" << static_cast<Real> (super::M_flag) << "): " << super::M_operator->flux( super::M_flag )  << "\n";
#endif

				setVariable( "f_flux", super::M_operator->flux( super::M_flag ) );

				break;

			case super::f_pressure :

#ifdef DEBUG
				Debug( 5027 ) << "                                               f_pressure(" << static_cast<Real> (super::M_flag) << "): " << super::M_operator->pressure( super::M_flag )  << "\n";
#endif

				setVariable( "f_pressure", super::M_operator->pressure( super::M_flag ) );

				break;
		}
}

} // Namespace LifeV

#endif /* __BCInterfaceOseenFunction_H */
