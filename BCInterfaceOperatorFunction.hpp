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





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/life.hpp>
#include <life/lifefem/bcFunction.hpp>

#include <string>

#include <lifemc/lifefem/BCInterfaceData.hpp>
#include <lifemc/lifefem/BCInterfaceFunction.hpp>





// ===================================================
//! Namespaces & Enums
// ===================================================
namespace LifeV {





/*!
 * \class BCInterfaceOperatorFunction
 * \brief LifeV bcFunction wrapper for BCInterface (with operators).
 *
 *  @author Cristiano Malossi
 *  @see
 *
 *  This class is an interface between BCInterface, SpiritParser and a general
 *  LifeV operator (such as Oseen or FSIOperator). It allows to construct LifeV
 *  functions type for boundary conditions, using a functions string loaded from
 *  a GetPot file in which are present some operator parameters.
 *
 *  The class can be used in two ways:
 *
 *  1) hereditating it with and implementing createAccessList(), addOperatorVariables();
 *  2) manually setting the variables by using the setVariable() function.
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
class BCInterfaceOperatorFunction : public virtual BCInterfaceFunction<Operator>
//     :
//     public LifeV::Application
{
public:

	typedef BCInterfaceFunction<Operator>			super;

	/** @name Constructors & Destructor
     */
    //@{

    //! Constructor
	BCInterfaceOperatorFunction();

	//! Constructor
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	BCInterfaceOperatorFunction( const BCInterfaceData<Operator>& data );

	//! Copy constructor
	/*!
	 * \param function			- BCInterfaceOperatorFunction
	 */
	BCInterfaceOperatorFunction( const BCInterfaceOperatorFunction& function );

    //! Destructor
    ~BCInterfaceOperatorFunction() {}

    //@}


	/** @name Methods
     */
    //@{

	//! Operator =
	/*!
	 * \param function			- BCInterfaceOperatorFunction
	 */
	virtual BCInterfaceOperatorFunction& operator=( const BCInterfaceOperatorFunction& function );

	//! Set data
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	virtual void setData( const BCInterfaceData<Operator>& data );

	//! Compare function
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	virtual bool compare( const BCInterfaceData<Operator>& data );

	//! Set variable function
	/*!
	 * \param name				- name of the variable
	 * \param value				- value of the variable
	 */
	inline void setVariable( const std::string& name, const Real& value );

	//@}

protected:

	/** @name Protected functions
	*/
	//@{

	virtual inline void createAccessList( void );

	virtual inline void addOperatorVariables( const Real& /*t*/ ) {}

	//@}

	//List of all available operators
	enum operatorList
	{
		f_area,
		f_flux,
		f_pressure,
		s_density,
		s_poisson,
		s_thickness,
		s_young,
	};

	boost::shared_ptr<Operator>							M_operator;
	BCFlag												M_flag;
	std::set<operatorList>								M_list;
	std::map<std::string, operatorList>					M_mapList;
	Real												M_oldTime;
};

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
	M_mapList						( ),
	M_oldTime						( -1.0 ) // Negative time
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
	M_mapList						( ),
	M_oldTime						( -1.0 ) // Negative time
{

#ifdef DEBUG
	Debug( 5023 ) << "BCInterfaceOperatorFunction::BCInterfaceOperatorFunction( data )" << "\n";
#endif

	this->setData( data );
}



template <class Operator>
BCInterfaceOperatorFunction<Operator>::BCInterfaceOperatorFunction( const BCInterfaceOperatorFunction& function ) :
	BCInterfaceFunction<Operator>	( function ),
	M_operator						( function.M_operator ),
	M_flag							( function.M_flag ),
	M_list							( function.M_list ),
	M_mapList						( function.M_mapList ),
	M_oldTime						( function.M_oldTime )
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
    	M_oldTime		= function.M_oldTime;
    }

	return *this;
}



template <class Operator>
void
BCInterfaceOperatorFunction<Operator>::setData( const BCInterfaceData<Operator>& data )
{

#ifdef DEBUG
	Debug( 5023 ) << "BCInterfaceOperatorFunction::setData" << "\n";
#endif

	M_operator	= data.get_operator();
	M_flag		= data.get_flag();

	super::setData( data );

	this->createAccessList();
}



template <class Operator>
bool
BCInterfaceOperatorFunction<Operator>::compare( const BCInterfaceData<Operator>& data )
{
	return	super::M_baseString.compare( data.get_baseString() ) == 0	&&
			super::M_comV == data.get_comV() && M_flag == data.get_flag();
}



template <class Operator>
inline void
BCInterfaceOperatorFunction<Operator>::setVariable( const std::string& name, const Real& value )
{
	super::M_parser->setVariable( name, value );
}





// ===================================================
//! Protected functions
// ===================================================
template <class Operator>
inline void
BCInterfaceOperatorFunction<Operator>::createAccessList( void )
{

#ifdef DEBUG
	Debug( 5023 ) << "BCInterfaceOperatorFunction::createAccessList" << "\n";
#endif

	//Create list
	M_list.clear();
	for ( typename std::map<std::string, operatorList>::iterator j = M_mapList.begin() ; j != M_mapList.end() ; ++j )
		if ( boost::find_first( super::M_baseString, j->first ) )
			M_list.insert( j->second );
}

} // Namespace LifeV

#endif /* __BCInterfaceOperatorFunction_H */
