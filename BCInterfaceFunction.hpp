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

#include <lifemc/lifecore/SpiritParser.hpp>





// ===================================================
//! Namespaces
// ===================================================
using namespace LifeV;





/*!
 * \class BCInterfaceFunction
 * \brief LifeV bcFunction wrapper for BCInterface.
 *
 *  @author Cristiano Malossi
 *  @see
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

    //! Constructor
	BCInterfaceFunction( const std::string& baseString, const std::string& stringSeparator=";" );

    //! Destructor
    ~BCInterfaceFunction() {}

    //@}



    /** @name Get functions
     */
    //@{

	BCFunctionBase& getBase() { return M_base; }

    //@}

private:

	// ===================================================
	//! Private variables
	// ===================================================

	std::map<UInt,UInt>									M_mapID;
	BCFunctionBase 										M_base;
	std::vector< boost::shared_ptr<SpiritParser> >		M_parserVector;


	// ===================================================
	//! Private functions
	// ===================================================

	/** @name Private functions
	*/
	//@{

    //! buildFunctionBase
    void buildFunctionBase( void );

    //! getFunction
    function_type getFunction( void );

    //! Function
    Real Function( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id );

    //@}
};

#endif /* __BCInterfaceFunction_H */
