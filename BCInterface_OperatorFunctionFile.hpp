//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief BCInterface_OperatorFunctionFile
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 26-08-2009
 */

#ifndef BCInterface_OperatorFunctionFile_H
#define BCInterface_OperatorFunctionFile_H 1

#include <lifemc/lifesolver/BCInterface_FunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface_OperatorFunction.hpp>

namespace LifeV {

//! BCInterface_OperatorFunctionFile - LifeV bcFunction wrapper for BCInterface (with Operators)
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface, SpiritParser and and a general
 *  LifeV operator (such as Oseen or FSIOperator). It allows to construct LifeV
 *  functions type for boundary conditions, using a GetPot file containing a function string and a
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
 *  The function string can contain Operator parameters.
 */
template< class Operator >
class BCInterface_OperatorFunctionFile: public virtual BCInterface_FunctionFile< Operator > ,
                                        public virtual BCInterface_OperatorFunction< Operator >
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface_Function< Operator >                        super0;
    typedef BCInterface_FunctionFile< Operator >                    super1;
    typedef BCInterface_OperatorFunction< Operator >                super2;

    typedef typename super0::Data_Type                              Data_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface_OperatorFunctionFile();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface_OperatorFunctionFile( const Data_Type& data );

    //! Copy constructor
    /*!
     * @param function BCInterface_OperatorFunctionFile
     */
    BCInterface_OperatorFunctionFile( const BCInterface_OperatorFunctionFile& function );

    //! Destructor
    virtual ~BCInterface_OperatorFunctionFile() {}

    //@}


    //! @name Methods
    //@{

    //! Operator =
    /*!
     * @param function BCInterface_OperatorFunctionFile
     * @return reference to a copy of the class
     */
    virtual BCInterface_OperatorFunctionFile& operator=( const BCInterface_OperatorFunctionFile& function );

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void SetData( const Data_Type& data );

    //@}

};

//! Factory create function
template< typename Operator >
inline BCInterface_Function< Operator >* BCInterface_CreateOperatorFunctionFile()
{
    return new BCInterface_OperatorFunctionFile< Operator > ();
}

// ===================================================
// Constructors
// ===================================================
template< class Operator >
BCInterface_OperatorFunctionFile< Operator >::BCInterface_OperatorFunctionFile() :
    super0      (),
    super1      (),
    super2      ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface_OperatorFunctionFile::BCInterface_OperatorFunctionFile()" << "\n";
#endif

}

template< class Operator >
BCInterface_OperatorFunctionFile< Operator >::BCInterface_OperatorFunctionFile( const Data_Type& data ) :
    super0      (),
    super1      (),
    super2      ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface_OperatorFunctionFile::BCInterface_OperatorFunctionFile( data )" << "\n";
#endif

    this->SetData( data );
}

template< class Operator >
BCInterface_OperatorFunctionFile< Operator >::BCInterface_OperatorFunctionFile( const BCInterface_OperatorFunctionFile& function ) :
    super0      ( function ),
    super1      ( function ),
    super2      ( function )
{
}

// ===================================================
// Methods
// ===================================================
template< class Operator >
BCInterface_OperatorFunctionFile< Operator >&
BCInterface_OperatorFunctionFile< Operator >::operator=( const BCInterface_OperatorFunctionFile& function )
{
    if ( this != &function )
    {
        super1::operator=( function );
        super2::operator=( function );
    }

    return *this;
}

template< class Operator >
void
BCInterface_OperatorFunctionFile< Operator >::SetData( const Data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface_OperatorFunctionFile::setData" << "\n";
#endif
    super1::SetData( data );

    //super2::SetData( data ); Cannot call directly, because it call again BCInterface_Function::setup( data )
    super2::M_operator = data.GetOperator();
    super2::M_flag     = data.GetFlag();

    super2::CreateAccessList( data );
}

} // Namespace LifeV

#endif /* BCInterface_OperatorFunctionFile_H */
