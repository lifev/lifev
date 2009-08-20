/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-07-22

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
   \file BCInterfaceFSIFunctionFile.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-07-22
 */

#ifndef __BCInterfaceFSIFunctionFile_H
#define __BCInterfaceFSIFunctionFile_H 1





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>

#include <sstream>

#include <lifemc/lifefem/BCInterfaceFunctionFile.hpp>
#include <lifemc/lifefem/BCInterfaceFSIFunction.hpp>





// ===================================================
//! Namespaces & Enums
// ===================================================
namespace LifeV {





/*!
 * \class BCInterfaceFSIFunctionFile
 * \brief LifeV bcFunction wrapper for BCInterface (FSI problems).
 *
 *  @author Cristiano Malossi
 *  @see
 *
 *  This class is an interface between BCInterface, SpiritParser and FSIOperator. It allows to construct LifeV
 *  functions type for boundary conditions, using a GetPot file containing a function string and a
 *  table of discrete data (for example a Flux or a Pressure depending on time). Moreover the function string can
 *  contain FSI parameters such as fluid flux and pressure, or solid density and thickness.
 *
 */
class BCInterfaceFSIFunctionFile : 	public BCInterfaceFunctionFile,
									public BCInterfaceFSIFunction
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
	BCInterfaceFSIFunctionFile( const boost::shared_ptr<FSIOperator>& Oper );

    //! Constructor
	/*!
	 * \param data				- BC data loaded from GetPot file
	 * \param Oper				- FSIOperator
	 */
	BCInterfaceFSIFunctionFile( const BCInterfaceData& data, const boost::shared_ptr<FSIOperator>& Oper );

	//! Copy constructor
	/*!
	 * \param function			- BCInterfaceFSIFunctionFile
	 */
	BCInterfaceFSIFunctionFile( const BCInterfaceFSIFunctionFile& function );

    //! Destructor
    ~BCInterfaceFSIFunctionFile() {}

    //@}



    /** @name Methods
     */
    //@{

    //! Operator =
    /*!
     * \param function			- BCInterfaceFSIFunctionFile
     */
    BCInterfaceFSIFunctionFile& operator=( const BCInterfaceFSIFunctionFile& function );

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

};

} // Namespace LifeV

#endif /* __BCInterfaceFSIFunctionFile_H */
