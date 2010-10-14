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
 *  @date 10-05-2010
 */

#ifndef BCInterface1D_OperatorFunctionFile_H
#define BCInterface1D_OperatorFunctionFile_H 1

#include <lifemc/lifesolver/BCInterface1D_FunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface1D_OperatorFunction.hpp>

namespace LifeV {

//! BCInterface1D_OperatorFunctionFile - LifeV bcFunction wrapper for BCInterface1D (with Operators)
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface1D, SpiritParser and and a general
 *  LifeV operator (such as Oseen or FSIOperator). It allows to construct LifeV
 *  functions type for boundary conditions, using a GetPot file containing a function string and a
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
 *  The function string can contain Operator parameters.
 */
template< class Operator >
class BCInterface1D_OperatorFunctionFile: public virtual BCInterface1D_FunctionFile< Operator > ,
                                          public virtual BCInterface1D_OperatorFunction< Operator >
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface1D_Function< Operator >                      super0;
    typedef BCInterface1D_FunctionFile< Operator >                  super1;
    typedef BCInterface1D_OperatorFunction< Operator >              super2;

    typedef BCInterface1D_Data                                      Data_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface1D_OperatorFunctionFile();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface1D_OperatorFunctionFile( const Data_Type& data );

    //! Copy constructor
    /*!
     * @param function BCInterface1D_OperatorFunctionFile
     */
    BCInterface1D_OperatorFunctionFile( const BCInterface1D_OperatorFunctionFile& function );

    //! Destructor
    virtual ~BCInterface1D_OperatorFunctionFile() {}

    //@}


    //! @name Methods
    //@{

    //! Operator =
    /*!
     * @param function BCInterface1D_OperatorFunctionFile
     * @return reference to a copy of the class
     */
    virtual BCInterface1D_OperatorFunctionFile& operator=( const BCInterface1D_OperatorFunctionFile& function );

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void SetData( const Data_Type& data );

    //@}

};

//! Factory create function
template< typename Operator >
inline BCInterface1D_Function< Operator >* BCInterface1D_CreateOperatorFunctionFile()
{
    return new BCInterface1D_OperatorFunctionFile< Operator > ();
}

// ===================================================
// Constructors
// ===================================================
template< class Operator >
BCInterface1D_OperatorFunctionFile< Operator >::BCInterface1D_OperatorFunctionFile() :
    super0      (),
    super1      (),
    super2      ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface1D_OperatorFunctionFile::BCInterface1D_OperatorFunctionFile()" << "\n";
#endif

}

template< class Operator >
BCInterface1D_OperatorFunctionFile< Operator >::BCInterface1D_OperatorFunctionFile( const Data_Type& data ) :
    super0      (),
    super1      (),
    super2      ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface1D_OperatorFunctionFile::BCInterface1D_OperatorFunctionFile( data )" << "\n";
#endif

    this->SetData( data );
}

template< class Operator >
BCInterface1D_OperatorFunctionFile< Operator >::BCInterface1D_OperatorFunctionFile( const BCInterface1D_OperatorFunctionFile& function ) :
    super0      ( function ),
    super1      ( function ),
    super2      ( function )
{
}

// ===================================================
// Methods
// ===================================================
template< class Operator >
BCInterface1D_OperatorFunctionFile< Operator >&
BCInterface1D_OperatorFunctionFile< Operator >::operator=( const BCInterface1D_OperatorFunctionFile& function )
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
BCInterface1D_OperatorFunctionFile< Operator >::SetData( const Data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface1D_OperatorFunctionFile::setData" << "\n";
#endif
    super1::SetData( data );

    //super2::SetData( data ); Cannot call directly, because it call again BCInterface1D1D_Function::setup( data )
    super2::M_side     = data.GetSide();

    super2::CreateAccessList( data );
}

} // Namespace LifeV

#endif /* BCInterface1D_OperatorFunctionFile_H */
