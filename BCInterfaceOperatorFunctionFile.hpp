/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-08-26

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
   \file BCInterfaceOperatorFunctionFile.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-08-26
 */

#ifndef __BCInterfaceOperatorFunctionFile_H
#define __BCInterfaceOperatorFunctionFile_H 1





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/life.hpp>

#include <lifemc/lifefem/BCInterfaceFunctionFile.hpp>
#include <lifemc/lifefem/BCInterfaceOperatorFunction.hpp>





// ===================================================
//! Namespaces & Enums
// ===================================================
namespace LifeV {





/*!
 * \class BCInterfaceOperatorFunctionFile
 * \brief LifeV bcFunction wrapper for BCInterface (with Operators).
 *
 *  @author Cristiano Malossi
 *  @see
 *
 *  This class is an interface between BCInterface, SpiritParser and and a general
 *  LifeV operator (such as Oseen or FSIOperator). It allows to construct LifeV
 *  functions type for boundary conditions, using a GetPot file containing a function string and a
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
 *  The function string can contain Operator parameters.
 *
 */
template <class Operator>
class BCInterfaceOperatorFunctionFile : public virtual BCInterfaceFunctionFile<Operator>,
										public virtual BCInterfaceOperatorFunction<Operator>
//     :
//     public LifeV::Application
{
public:

	typedef BCInterfaceFunctionFile<Operator>		super1;
	typedef BCInterfaceOperatorFunction<Operator>	super2;

	/** @name Constructors & Destructor
     */
    //@{

    //! Constructor
	BCInterfaceOperatorFunctionFile();

    //! Constructor
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	BCInterfaceOperatorFunctionFile( const BCInterfaceData<Operator>& data );

	//! Copy constructor
	/*!
	 * \param function			- BCInterfaceOperatorFunctionFile
	 */
	BCInterfaceOperatorFunctionFile( const BCInterfaceOperatorFunctionFile& function );

    //! Destructor
    ~BCInterfaceOperatorFunctionFile() {}

    //@}



    /** @name Methods
     */
    //@{

    //! Operator =
    /*!
     * \param function			- BCInterfaceOperatorFunctionFile
     */
    virtual BCInterfaceOperatorFunctionFile& operator=( const BCInterfaceOperatorFunctionFile& function );

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

    //@}

};

//! Factory create function
template <typename Operator>
inline BCInterfaceFunction<Operator>* createOperatorFunctionFile()
{
	return new BCInterfaceOperatorFunctionFile<Operator>();
}



// ===================================================
//! Constructors
// ===================================================
template <class Operator>
BCInterfaceOperatorFunctionFile<Operator>::BCInterfaceOperatorFunctionFile( ) :
	BCInterfaceFunction<Operator>				( ),
	BCInterfaceFunctionFile<Operator>			( ),
	BCInterfaceOperatorFunction<Operator>		( )
{

#ifdef DEBUG
	Debug( 5024 ) << "BCInterfaceOperatorFunctionFile::BCInterfaceOperatorFunctionFile( void )" << "\n";
#endif

}



template <class Operator>
BCInterfaceOperatorFunctionFile<Operator>::BCInterfaceOperatorFunctionFile( const BCInterfaceData<Operator>& data ) :
	BCInterfaceFunction<Operator>				( ),
	BCInterfaceFunctionFile<Operator>			( ),
	BCInterfaceOperatorFunction<Operator>		( )
{

#ifdef DEBUG
	Debug( 5024 ) << "BCInterfaceOperatorFunctionFile::BCInterfaceOperatorFunctionFile( data )" << "\n";
#endif

	this->setData( data );
}



template <class Operator>
BCInterfaceOperatorFunctionFile<Operator>::BCInterfaceOperatorFunctionFile( const BCInterfaceOperatorFunctionFile& function ) :
	BCInterfaceFunction<Operator>				( function ),
	BCInterfaceFunctionFile<Operator>			( function ),
	BCInterfaceOperatorFunction<Operator>		( function )
{
}





// ===================================================
//! Methods
// ===================================================
template <class Operator>
BCInterfaceOperatorFunctionFile<Operator>&
BCInterfaceOperatorFunctionFile<Operator>::operator=( const BCInterfaceOperatorFunctionFile& function )
{
	if ( this != &function )
	{
		super1::operator=( function );
		super2::operator=( function );
	}

	return *this;
}



template <class Operator>
void
BCInterfaceOperatorFunctionFile<Operator>::setData( const BCInterfaceData<Operator>& data )
{

#ifdef DEBUG
	Debug( 5024 ) << "BCInterfaceOperatorFunctionFile::setData" << "\n";
#endif

	super1::setData( data );
	super2::setData( data );
}



template <class Operator>
bool
BCInterfaceOperatorFunctionFile<Operator>::compare( const BCInterfaceData<Operator>& data )
{
	return super1::M_fileName.compare( data.get_baseString() ) == 0 && super1::M_comV == data.get_comV() && super2::M_flag == data.get_flag();
}

} // Namespace LifeV

#endif /* __BCInterfaceOperatorFunctionFile_H */
