/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-08-24

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
   \file BCInterfaceOperatorFunction.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-08-24
 */

#ifndef __BCInterfaceOperatorFunction_H
#define __BCInterfaceOperatorFunction_H 1

#include <life/lifecore/life.hpp>
#include <life/lifefem/bcFunction.hpp>

#include <string>

#include <life/lifesolver/FSIOperator.hpp>
#include <life/lifesolver/Oseen.hpp>

#include <lifemc/lifefem/BCInterfaceData.hpp>
#include <lifemc/lifefem/BCInterfaceFunction.hpp>

namespace LifeV {


//! BCInterfaceOperatorFunction - LifeV bcFunction wrapper for BCInterface (with operators)
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface, SpiritParser and a general
 *  LifeV operator (such as Oseen or FSIOperator). It allows to construct LifeV
 *  functions type for boundary conditions, using a functions string loaded from
 *  a GetPot file in which are present some operator parameters.
 *
 *  The class can be used in two ways:
 *
 *  1) hereditating it and implementing the template specialization of createAccessList() and UpdateOperatorVariables();
 *  2) manually setting the variables by using the setVariable() function.
 *
 *	<b>AVAILABLE OPERATORS</b>
 *
 *	Available operators are:
 *
 *	f_area
 *	f_density
 *	f_flux
 *	f_pressure
 *	f_viscosity
 *	s_density
 *	s_poisson
 *	s_thickness
 *	s_young
 */
template <class Operator>
class BCInterfaceOperatorFunction : public virtual BCInterfaceFunction<Operator>
{
public:

	typedef BCInterfaceFunction<Operator>			super;

	//! @name Constructors & Destructor
    //@{

    //! Constructor
	BCInterfaceOperatorFunction();

	//! Constructor
	/*!
	 * \param data - BC data loaded from GetPot file
	 */
	BCInterfaceOperatorFunction( const BCInterfaceData<Operator>& data );

	//! Copy constructor
	/*!
	 * \param function - BCInterfaceOperatorFunction
	 */
	BCInterfaceOperatorFunction( const BCInterfaceOperatorFunction& function );

    //! Destructor
	virtual ~BCInterfaceOperatorFunction() {}

    //@}



	//! @name Methods
    //@{

	//! Operator =
	/*!
	 * \param function - BCInterfaceOperatorFunction
	 */
	virtual BCInterfaceOperatorFunction& operator=( const BCInterfaceOperatorFunction& function );

	//! Set data
	/*!
	 * \param data - BC data loaded from GetPot file
	 */
	virtual void SetData( const BCInterfaceData<Operator>& data );

	//! Compare function
	/*!
	 * \param data - BC data loaded from GetPot file
	 */
	virtual bool Compare( const BCInterfaceData<Operator>& data );

	//! Set variable function
	/*!
	 * \param name - name of the variable
	 * \param value - value of the variable
	 */
	inline void SetVariable( const std::string& name, const Real& value );

	//! Update operator variables
	inline void UpdateOperatorVariables( void ) {}

	//@}

protected:

	//! @name Protected functions
	//@{

	inline void CreateAccessList( void );

	//@}

	//List of all available operators
	enum operatorList
	{
		f_area,
		f_flux,
		f_density,
		f_pressure,
		f_viscosity,
		s_density,
		s_poisson,
		s_thickness,
		s_young,
	};

	boost::shared_ptr<Operator>							M_operator;
	BCFlag												M_flag;
	std::set<operatorList>								M_list;
	std::map<std::string, operatorList>					M_mapList;

private:

	inline void CreateFluidMap( void );
	inline void CreateSolidMap( void );
	inline void CreateList    ( void );

    inline void SwitchErrorMessage( const std::string& operatorType );

};

//! Factory create function
template <typename Operator>
inline BCInterfaceFunction<Operator>* createOperatorFunction()
{
	return new BCInterfaceOperatorFunction<Operator>();
}



// ===================================================
//! Constructors
// ===================================================
template <class Operator>
BCInterfaceOperatorFunction<Operator>::BCInterfaceOperatorFunction( ) :
	BCInterfaceFunction<Operator>	( ),
	M_operator						( ),
	M_flag							( ),
	M_list							( ),
	M_mapList						( )
{

#ifdef DEBUG
	Debug( 5023 ) << "BCInterfaceOperatorFunction::BCInterfaceOperatorFunction( void )" << "\n";
#endif

}

template <class Operator>
BCInterfaceOperatorFunction<Operator>::BCInterfaceOperatorFunction( const BCInterfaceData<Operator>& data ) :
	BCInterfaceFunction<Operator>	( ),
	M_operator						( ),
	M_flag							( ),
	M_list							( ),
	M_mapList						( )
{

#ifdef DEBUG
	Debug( 5023 ) << "BCInterfaceOperatorFunction::BCInterfaceOperatorFunction( data )" << "\n";
#endif

	this->SetData( data );
}

template <class Operator>
BCInterfaceOperatorFunction<Operator>::BCInterfaceOperatorFunction( const BCInterfaceOperatorFunction& function ) :
	BCInterfaceFunction<Operator>	( function ),
	M_operator						( function.M_operator ),
	M_flag							( function.M_flag ),
	M_list							( function.M_list ),
	M_mapList						( function.M_mapList )
{
}



// ===================================================
//! Methods
// ===================================================
template <class Operator>
BCInterfaceOperatorFunction<Operator>&
BCInterfaceOperatorFunction<Operator>::operator=( const BCInterfaceOperatorFunction& function )
{
    if ( this != &function )
    {
    	super::operator=( function );

    	M_operator		= function.M_operator;
    	M_flag			= function.M_flag;
    	M_list			= function.M_list;
    	M_mapList		= function.M_mapList;
    }

	return *this;
}

template <class Operator>
void
BCInterfaceOperatorFunction<Operator>::SetData( const BCInterfaceData<Operator>& data )
{

#ifdef DEBUG
	Debug( 5023 ) << "BCInterfaceOperatorFunction::setData" << "\n";
#endif

	M_operator	= data.GetOperator();
	M_flag		= data.GetFlag();

	super::SetData( data );

	CreateAccessList();
}

template <class Operator>
bool
BCInterfaceOperatorFunction<Operator>::Compare( const BCInterfaceData<Operator>& data )
{
	return	super::M_baseString.compare( data.GetBaseString() ) == 0	&&
			super::M_comV == data.GetComV() && M_flag == data.GetFlag();
}

template <class Operator>
inline void
BCInterfaceOperatorFunction<Operator>::SetVariable( const std::string& name, const Real& value )
{
	super::M_parser->setVariable( name, value );
}

template <>
inline void
BCInterfaceOperatorFunction<FSIOperator>::UpdateOperatorVariables( void)
{

#ifdef DEBUG
	Debug( 5023 ) << "BCInterfaceOperatorFunction<FSIOperator>::UpdateOperatorVariables  " << "\n";
#endif

	// Create/Update variables for FSI problem
	for ( std::set<operatorList>::iterator j = M_list.begin() ; j != M_list.end() ; ++j )
		switch ( *j )
		{
			// f_ -> FLUID
			case f_area :

#ifdef DEBUG
				Debug( 5023 ) << "                                                   f_area(" << static_cast<Real> (M_flag) << "): " << M_operator->fluid().area( M_flag )  << "\n";
#endif
				SetVariable( "f_area", M_operator->fluid().area( M_flag ) );

				break;

			case f_density :

#ifdef DEBUG
				Debug( 5023 ) << "                                                f_density(): " << M_operator->fluid().density()  << "\n";
#endif
				SetVariable( "f_density", M_operator->fluid().density() );

				break;

			case f_flux :

#ifdef DEBUG
				Debug( 5023 ) << "                                                   f_flux(" << static_cast<Real> (M_flag) << "): " << M_operator->fluid().flux( M_flag )  << "\n";
#endif

				SetVariable( "f_flux", M_operator->fluid().flux( M_flag ) );

				break;

			case f_pressure :

#ifdef DEBUG
				Debug( 5023 ) << "                                               f_pressure(" << static_cast<Real> (M_flag) << "): " << M_operator->fluid().pressure( M_flag )  << "\n";
#endif

				SetVariable( "f_pressure", M_operator->fluid().pressure( M_flag ) );

				break;

			case f_viscosity :

#ifdef DEBUG
				Debug( 5023 ) << "                                              f_viscosity(): " << M_operator->fluid().viscosity()  << "\n";
#endif
				SetVariable( "f_viscosity", M_operator->fluid().viscosity() );

				break;


			// s_ -> SOLID
			case s_density :

#ifdef DEBUG
				Debug( 5023 ) << "                                                   s_density: " << M_operator->solid().rho()  << "\n";
#endif

				SetVariable( "s_density", M_operator->solid().rho() );

				break;

			case s_poisson :

#ifdef DEBUG
				Debug( 5023 ) << "                                                   s_poisson: " << M_operator->solid().poisson()  << "\n";
#endif

				SetVariable( "s_poisson", M_operator->solid().poisson() );

				break;

			case s_thickness :

#ifdef DEBUG
				Debug( 5023 ) << "                                                 s_thickness: " << M_operator->solid().thickness()  << "\n";
#endif

				SetVariable( "s_thickness", M_operator->solid().thickness() );

				break;

			case s_young :

#ifdef DEBUG
				Debug( 5023 ) << "                                                     s_young: " << M_operator->solid().young()  << "\n";
#endif

				SetVariable( "s_young", M_operator->solid().young() );

				break;

			default :
				SwitchErrorMessage( "FSI" );
		}
}

template <>
inline void
BCInterfaceOperatorFunction< Oseen< RegionMesh3D< LinearTetra > > >::UpdateOperatorVariables( void )
{

#ifdef DEBUG
	Debug( 5023 ) << "BCInterfaceOperatorFunction<Oseen>::UpdateOperatorVariables  " << "\n";
#endif

	// Create/Update variables for FSI problem
	for ( std::set<operatorList>::iterator j = M_list.begin() ; j != M_list.end() ; ++j )
		switch ( *j )
		{
			// f_ -> FLUID
			case f_area :

#ifdef DEBUG
				Debug( 5023 ) << "                                                   f_area(" << static_cast<Real> (M_flag) << "): " << M_operator->area( M_flag )  << "\n";
#endif
				SetVariable( "f_area", M_operator->area( M_flag ) );

				break;

			case f_density :

#ifdef DEBUG
				Debug( 5023 ) << "                                                f_density(): " << M_operator->density()  << "\n";
#endif
				SetVariable( "f_density", M_operator->density() );

				break;

			case f_flux :

#ifdef DEBUG
				Debug( 5023 ) << "                                                   f_flux(" << static_cast<Real> (M_flag) << "): " << M_operator->flux( M_flag )  << "\n";
#endif

				SetVariable( "f_flux", M_operator->flux( M_flag ) );

				break;

			case f_pressure :

#ifdef DEBUG
				Debug( 5023 ) << "                                               f_pressure(" << static_cast<Real> (M_flag) << "): " << M_operator->pressure( M_flag )  << "\n";
#endif

				SetVariable( "f_pressure", M_operator->pressure( M_flag ) );

				break;

			case f_viscosity :

#ifdef DEBUG
				Debug( 5023 ) << "                                              f_viscosity(): " << M_operator->viscosity()  << "\n";
#endif
				SetVariable( "f_viscosity", M_operator->viscosity() );

				break;

			default :

				SwitchErrorMessage( "OSEEN" );
		}
}



// ===================================================
//! Protected functions
// ===================================================
template <>
inline void
BCInterfaceOperatorFunction<FSIOperator>::CreateAccessList( void )
{

#ifdef DEBUG
	Debug( 5023 ) << "BCInterfaceOperatorFunction<FSIOperator>::createAccessList" << "\n";
#endif

	CreateFluidMap();
	CreateSolidMap();
	CreateList();
}

template <>
inline void
BCInterfaceOperatorFunction< Oseen< RegionMesh3D< LinearTetra > > >::CreateAccessList( void )
{

#ifdef DEBUG
	Debug( 5023 ) << "BCInterfaceOperatorFunction<Oseen>::createAccessList" << "\n";
#endif

	CreateFluidMap();
	CreateList();
}



// ===================================================
//! Private functions
// ===================================================
template <class Operator>
inline void
BCInterfaceOperatorFunction<Operator>::CreateFluidMap( void )
{
	M_mapList["f_area"]		= f_area;
	M_mapList["f_density"]	= f_density;
	M_mapList["f_flux"]		= f_flux;
	M_mapList["f_pressure"]	= f_pressure;
	M_mapList["f_viscosity"]= f_viscosity;
}

template <class Operator>
inline void
BCInterfaceOperatorFunction<Operator>::CreateSolidMap( void )
{
	M_mapList["s_density"]	= s_density;
	M_mapList["s_poisson"]	= s_poisson;
	M_mapList["s_thickness"]= s_thickness;
	M_mapList["s_young"]	= s_young;
}

template <class Operator>
inline void
BCInterfaceOperatorFunction<Operator>::CreateList( void )
{
	M_list.clear();
	for ( typename std::map<std::string, operatorList>::iterator j = M_mapList.begin() ; j != M_mapList.end() ; ++j )
		if ( boost::find_first( super::M_baseString, j->first ) )
			M_list.insert( j->second );
}

template <class Operator>
inline void
BCInterfaceOperatorFunction<Operator>::SwitchErrorMessage( const std::string& operatorType )
{
	std::cout	<< "ERROR: Invalid variable type for " << operatorType << "OperatorFunction" << std::endl;
}

} // Namespace LifeV

#endif /* __BCInterfaceOperatorFunction_H */
