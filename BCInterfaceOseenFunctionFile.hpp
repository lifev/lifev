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
   \file BCInterfaceOseenFunctionFile.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-08-27
 */

#ifndef __BCInterfaceOseenFunctionFile_H
#define __BCInterfaceOseenFunctionFile_H 1





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/life.hpp>

#include <lifemc/lifefem/BCInterfaceOperatorFunctionFile.hpp>
#include <lifemc/lifefem/BCInterfaceOseenFunction.hpp>





// ===================================================
//! Namespaces & Enums
// ===================================================
namespace LifeV {





/*!
 * \class BCInterfaceOseenFunctionFile
 * \brief LifeV bcFunction wrapper for BCInterface (Oseen problems).
 *
 *  @author Cristiano Malossi
 *  @see
 *
 *  This class is a specialization of BCInterfaceOperatorFunctionFile class for Oseen problems.
 *
 */
template <class Operator>
class BCInterfaceOseenFunctionFile : 	public BCInterfaceOperatorFunctionFile<Operator>,
										public BCInterfaceOseenFunction<Operator>
//     :
//     public LifeV::Application
{
public:

	/** @name Constructors & Destructor
     */
    //@{

    //! Constructor
	BCInterfaceOseenFunctionFile();

    //! Constructor
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	BCInterfaceOseenFunctionFile( const BCInterfaceData<Operator>& data );

	//! Copy constructor
	/*!
	 * \param function			- BCInterfaceOseenFunctionFile
	 */
	BCInterfaceOseenFunctionFile( const BCInterfaceOseenFunctionFile& function );

    //! Destructor
    ~BCInterfaceOseenFunctionFile() {}

    //@}

};

//! Factory create function
template <> template <typename Mesh, typename SolverType>
inline BCInterfaceFunction< Oseen<Mesh, SolverType> >* createOperatorFunction< Oseen<Mesh, SolverType> >()
{
	return new BCInterfaceOseenFunctionFile< Oseen<Mesh, SolverType> >();
}



// ===================================================
//! Constructors
// ===================================================
template <class Operator>
BCInterfaceOseenFunctionFile<Operator>::BCInterfaceOseenFunctionFile( ) :
	BCInterfaceFunction<Operator>				( ),
	BCInterfaceFunctionFile<Operator>			( ),
	BCInterfaceOperatorFunction<Operator>		( ),
	BCInterfaceOperatorFunctionFile<Operator>	( ),
	BCInterfaceOseenFunction<Operator>			( )
{

#ifdef DEBUG
	Debug( 5026 ) << "BCInterfaceOseenFunctionFile::BCInterfaceOseenFunctionFile( void )" << "\n";
#endif

}



template <class Operator>
BCInterfaceOseenFunctionFile<Operator>::BCInterfaceOseenFunctionFile( const BCInterfaceData<Operator>& data ) :
	BCInterfaceFunction<Operator>				( ),
	BCInterfaceFunctionFile<Operator>			( ),
	BCInterfaceOperatorFunction<Operator>		( ),
	BCInterfaceOperatorFunctionFile<Operator>	( ),
	BCInterfaceOseenFunction<Operator>			( )
{

#ifdef DEBUG
	Debug( 5026 ) << "BCInterfaceOseenFunctionFile::BCInterfaceOseenFunctionFile( data )" << "\n";
#endif

	this->setData( data );
}



template <class Operator>
BCInterfaceOseenFunctionFile<Operator>::BCInterfaceOseenFunctionFile( const BCInterfaceOseenFunctionFile& function ) :
	BCInterfaceFunction<Operator>				( function ),
	BCInterfaceFunctionFile<Operator>			( function ),
	BCInterfaceOperatorFunction<Operator>		( function ),
	BCInterfaceOperatorFunctionFile<Operator>	( function ),
	BCInterfaceOseenFunction<Operator>			( function )
{
}

} // Namespace LifeV

#endif /* __BCInterfaceOseenFunctionFile_H */
