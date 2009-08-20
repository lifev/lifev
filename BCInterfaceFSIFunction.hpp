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
#include <lifemc/lifefem/BCInterfaceFunction.hpp>





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
 *  This class is an interface between BCInterface, SpiritParser and FSIOperator. It allows to construct LifeV
 *  functions type for boundary conditions, using a functions string loaded from a GetPot file in which are present
 *  FSI parameters such as fluid flux and pressure, or solid density and thickness.
 *
 *	<b>FSI AVAILABLE OPERATORS</b>
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
class BCInterfaceFSIFunction : public virtual BCInterfaceFunction
//     :
//     public LifeV::Application
{
public:

	// ===================================================
	//! Public functions
	// ===================================================

	/** @name Constructors & Destructor
     */
    //@{

    //! Constructor
	/*!
	 * \param Oper				- FSIOperator
	 */
	BCInterfaceFSIFunction( const boost::shared_ptr<FSIOperator>& Oper );

    //! Constructor
	/*!
	 * \param data				- BC data loaded from GetPot file
	 * \param Oper				- FSIOperator
	 */
	BCInterfaceFSIFunction( const BCInterfaceData& data, const boost::shared_ptr<FSIOperator>& Oper );

	//! Copy constructor
	/*!
	 * \param function			- BCInterfaceFSIFunction
	 */
	BCInterfaceFSIFunction( const BCInterfaceFSIFunction& function );

    //! Destructor
    ~BCInterfaceFSIFunction() {}

    //@}


	/** @name Methods
     */
    //@{

	//! Operator =
	/*!
	 * \param function			- BCInterfaceFSIFunction
	 */
	BCInterfaceFSIFunction& operator=( const BCInterfaceFSIFunction& function );

	//! Set data
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	void setData( const BCInterfaceData& data );

	//! Compare function
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	bool compare( const BCInterfaceData& data );

	//@}

protected:

	enum FSIList
	{
		f_area,
		f_flux,
		f_pressure,
		s_density,
		s_poisson,
		s_thickness,
		s_young,
	};

	// ===================================================
	//! Private variables
	// ===================================================

	boost::shared_ptr<FSIOperator>						M_FSIOperator;
	BCFlag												M_flag;
	std::map<std::string,FSIList>						M_mapFSIList;
	std::set<FSIList>									M_FSIList;
	Real												M_oldTime;



	// ===================================================
	//! Private functions
	// ===================================================

	/** @name Private functions
	*/
	//@{

	inline void createFSIAccessList( void );

	void addFSIVariables( const Real& t );

	//@}
};

} // Namespace LifeV

#endif /* __BCInterfaceFSIFunction_H */
