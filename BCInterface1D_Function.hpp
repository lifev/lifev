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
 *  @brief BCInterface1D_Function
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 10-05-2010
 */

#ifndef BCInterface1D_Function_H
#define BCInterface1D_Function_H 1

#include <lifemc/lifesolver/BCInterface1D_Definitions.hpp>
#include <lifemc/lifesolver/BCInterface1D_Data.hpp>

#include <lifemc/lifecore/Parser.hpp>

namespace LifeV {

//! BCInterface1D_Function - LifeV bcFunction wrapper for BCInterface
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface and SpiritParser. It allows to construct LifeV
 *  functions type for boundary conditions, using a functions string loaded from a GetPot file.
 *
 *  <b>DETAILS:</b>
 *
 *  The constructor of the class takes a string contains the GetPot file function. By default the stringSeparator
 *  is set to semicolon ";".
 *
 *  The function string has to be in this form:
 *
 *  function = q
 *
 *  where q(t).
 *  Here there is an example:
 *
 *  function = '2*sin(2*pi*t)'
 *
 *  To set a constant for complicate expression it is possible to add them before the expression
 *  using a semicolon ";":
 *
 *  function = 'a=5.67436; (a*sin(2*pi*t))'
 */
template< typename Operator >
class BCInterface1D_Function
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface1D_Data                                                    Data_Type;
    typedef OneDimensionalModel_BCFunction                                        BCFunction_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    BCInterface1D_Function();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface1D_Function( const Data_Type& data );

    //! Copy constructor
    /*!
     * @param function BCInterface1D_Function
     */
    BCInterface1D_Function( const BCInterface1D_Function& function );

    //! Destructor
    virtual ~BCInterface1D_Function() {}

    //@}


    //! @name Methods
    //@{

    //! Operator =
    /*!
     * @param function BCInterface1D_Function
     * @return reference to a copy of the class
     */
    virtual BCInterface1D_Function& operator=( const BCInterface1D_Function& function );

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void SetData( const Data_Type& data );

    //@}


    //! @name Get Methods
    //@{

    //! Get the base of the boundary condition
    BCFunction_Type& GetBase();

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! dataInterpolation
    virtual inline void DataInterpolation() {}

    //@}

    boost::shared_ptr< Parser >       M_parser;

private:

    //! @name Private Methods
    //@{

    //! SetFunction
    inline void SetFunction();

    //! Function
    Real Function( const Real& t );

    //@}

    BCFunction_Type                   M_base;

};

//! Factory create function
template< typename Operator >
inline BCInterface1D_Function< Operator >* BCInterface1D_CreateFunction()
{
    return new BCInterface1D_Function< Operator > ();
}

// ===================================================
// Constructor
// ===================================================
template< typename Operator >
BCInterface1D_Function< Operator >::BCInterface1D_Function() :
    M_parser    (),
    M_base      ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface1D_Function::BCInterface1D_Function()" << "\n";
#endif

}

template< typename Operator >
BCInterface1D_Function< Operator >::BCInterface1D_Function( const Data_Type& data ) :
    M_parser    (),
    M_base      ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface1D_Function::BCInterface1D_Function( data )" << "\n";
#endif

    this->SetData( data );
}

template< typename Operator >
BCInterface1D_Function< Operator >::BCInterface1D_Function( const BCInterface1D_Function& function ) :
    M_parser    ( function.M_parser ),
    M_base      ( function.M_base )
{
}

// ===================================================
// Methods
// ===================================================
template< typename Operator >
BCInterface1D_Function< Operator >&
BCInterface1D_Function< Operator >::operator=( const BCInterface1D_Function& function )
{
    if ( this != &function )
    {
        M_parser     = function.M_parser;
        M_base       = function.M_base;
    }

    return *this;
}

template< typename Operator >
void
BCInterface1D_Function< Operator >::SetData( const Data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5022 ) << "BCInterface1D_Function::setData" << "\n";
#endif

    if ( M_parser )
        M_parser->SetString( data.GetBaseString() );
    else
        M_parser.reset( new Parser( data.GetBaseString() ) );

    SetFunction();
}

// ===================================================
// Get Methods
// ===================================================
template< typename Operator >
typename BCInterface1D_Function< Operator >::BCFunction_Type&
BCInterface1D_Function< Operator >::GetBase()
{
    return M_base;
}

// ===================================================
// Private Methods
// ===================================================
template< typename Operator >
inline void
BCInterface1D_Function< Operator >::SetFunction()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface1D_Function::setFunction\n";
#endif

    M_base.setFunction( boost::bind( &BCInterface1D_Function::Function, this, _1 ) );
}

template< typename Operator >
Real
BCInterface1D_Function< Operator >::Function( const Real& t )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface1D_Function::Function: " << "\n";
    Debug( 5021 ) << "                                                           t: " << t << "\n";
#endif

    M_parser->SetVariable( "t", t );

    this->DataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "                                                evaluate(" << 1 << ") : " << M_parser->Evaluate( 1 ) << "\n";
#endif

    return M_parser->Evaluate( 1 );
}

} // Namespace LifeV

#endif /* BCInterface1D_Function_H */
