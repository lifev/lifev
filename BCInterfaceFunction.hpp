/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-04-06

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
   \file BCInterfaceFunction.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-04-06
 */

#ifndef __BCInterfaceFunction_H
#define __BCInterfaceFunction_H 1





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/life.hpp>
#include <life/lifefem/bcFunction.hpp>

#include <string>

#include <lifemc/lifefem/BCInterfaceData.hpp>
#include <lifemc/lifecore/SpiritParser.hpp>





// ===================================================
//! Namespaces & Enums
// ===================================================
namespace LifeV {





/*!
 * \class BCInterfaceFunction
 * \brief LifeV bcFunction wrapper for BCInterface.
 *
 *  @author Cristiano Malossi
 *  @see
 *
 *  This class is an interface between BCInterface and SpiritParser. It allows to construct LifeV
 *  functions type for boundary conditions, using a functions string loaded from a GetPot file.
 *
 *
 *
 *  <b>DETAILS:</b>
 *
 *  The constructor of the class takes a string contains the GetPot file function. By default the stringSeparator
 *  is set to semicolon ";".
 *
 *  The function string has to be in this form:
 *
 *  function = '(u, v, w)'
 *
 *  where u(x,y,z,t), v(x,y,z,t), w(x,y,z,t).
 *  Here there is an example:
 *
 *  function = '(x^2 + y^2, 0, 2*sin(2*pi*t))'
 *
 *  To set a constant for complicate expression it is possible to add them before the expression
 *  using a semicolon ";":
 *
 *  function = 'a=5.67436; (x^2+y^2,0,c*sin(2*pi*C*t)^C)'
 *
 *  NOTE:
 *  In the boundary condition file, if you have three component with the same expression
 *  (the same function) you can both write:
 *
 *  function = '(0, 0, 0)'
 *
 *  and
 *
 *  function = 0
 *
 *  The only difference is that the second kind of instruction is more efficient during execution.
 *
 */
class BCInterfaceFunction
//     :
//     public LifeV::Application
{
public:

	// ===================================================
	//! Typedef
	// ===================================================

	typedef boost::function<Real ( 	Real const& t,
									Real const& x,
									Real const& y,
									Real const& z,
									ID 	 const& id	)> 	function_type;



	// ===================================================
	//! Public functions
	// ===================================================

	/** @name Constructors & Destructor
     */
    //@{

	//! Empty Constructor
	BCInterfaceFunction( void );

	//! Constructor
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	BCInterfaceFunction( const BCInterfaceData& data );

	//! Copy constructor
	/*!
	 * \param function			- BCInterfaceFunction
	 */
	BCInterfaceFunction( const BCInterfaceFunction& function );

    //! Destructor
    ~BCInterfaceFunction() {}

    //@}



    /** @name Methods
     */
    //@{

    //! Operator =
    /*!
     * \param function			- BCInterfaceFunction
     */
    BCInterfaceFunction& operator=( const BCInterfaceFunction& function );

    //! Set data
    /*!
	 * \param data				- BC data loaded from GetPot file
	 */
    virtual void setData( const BCInterfaceData& data );

	//! Compare function
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
    virtual bool compare( const BCInterfaceData& data );

    //@}



    /** @name Get functions
     */
    //@{

	BCFunctionBase& getBase() { return M_base; }

    //@}

protected:

	// ===================================================
	//! Private variables
	// ===================================================

	std::string											M_baseString;
	BCComV												M_comV;
	BCFunctionBase 										M_base;
	boost::shared_ptr<SpiritParser>						M_parser;
	std::map<ID, ID>									M_mapID;



	// ===================================================
	//! Private functions
	// ===================================================

	/** @name Private functions
	 */
	//@{

	//! dataInterpolation
	virtual inline void dataInterpolation( void ) {}

	//! addFSIVariables
	virtual inline void addFSIVariables( const Real& /*t*/ ) {}

    //! SetFunction
    void setFunction( void );

    //! Function
    Real Function( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*id*/ );

    //! FunctionID
    Real FunctionID( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id );

    //@}
};

} // Namespace LifeV

#endif /* __BCInterfaceFunction_H */
