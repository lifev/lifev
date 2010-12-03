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
 *  @brief BCInterface_Function
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 06-04-2009
 */

#ifndef BCInterface_Function_H
#define BCInterface_Function_H 1

#include <lifemc/lifesolver/BCInterface_Definitions.hpp>
#include <lifemc/lifesolver/BCInterface_Data.hpp>

#include <lifemc/lifecore/Parser.hpp>

namespace LifeV
{

//! BCInterface_Function - LifeV bcFunction wrapper for BCInterface
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
 *  function = '[u, v, w]'
 *
 *  where u(x,y,z,t), v(x,y,z,t), w(x,y,z,t).
 *  Here there is an example:
 *
 *  function = '[x^2 + y^2, 0, 2*sin(2*pi*t)]'
 *
 *  To set a constant for complicate expression it is possible to add them before the expression
 *  using a semicolon ";":
 *
 *  function = 'a=5.67436; [x^2+y^2,0,a*sin(2*pi*t)^a]'
 *
 *  NOTE:
 *  In the boundary condition file, if you have three component with the same expression
 *  (the same function) you can both write:
 *
 *  function = '[0, 0, 0]'
 *
 *  and
 *
 *  function = 0
 *
 *  The only difference is that the second kind of instruction is more efficient during execution.
 *
 */
template< typename Operator >
class BCInterface_Function
{
public:

    //! @name Type definitions
    //@{

    typedef BCFunctionBase                                                        BCFunction_Type;
    typedef BCInterface_Data                                                      Data_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    BCInterface_Function();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface_Function( const Data_Type& data );

    //! Copy constructor
    /*!
     * @param function BCInterface_Function
     */
    BCInterface_Function( const BCInterface_Function& function );

    //! Destructor
    virtual ~BCInterface_Function() {}

    //@}


    //! @name Methods
    //@{

    //! Operator =
    /*!
     * @param function BCInterface_Function
     * @return reference to a copy of the class
     */
    virtual BCInterface_Function& operator=( const BCInterface_Function& function );

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

    //! Function
    Real Function( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*id*/);

    //! FunctionID
    Real FunctionID( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id );

    //@}

    BCFunction_Type                  M_base;
    std::map< ID, ID >               M_mapID;

};

//! Factory create function
template< typename Operator >
inline BCInterface_Function< Operator >* BCInterface_CreateFunction()
{
    return new BCInterface_Function< Operator > ();
}

// ===================================================
// Constructor
// ===================================================
template< typename Operator >
BCInterface_Function< Operator >::BCInterface_Function() :
        M_parser    (),
        M_base      (),
        M_mapID     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface_Function::BCInterface_Function()" << "\n";
#endif

}

template< typename Operator >
BCInterface_Function< Operator >::BCInterface_Function( const Data_Type& data ) :
        M_parser    (),
        M_base      (),
        M_mapID     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface_Function::BCInterface_Function( data )" << "\n";
#endif

    this->SetData( data );
}

template< typename Operator >
BCInterface_Function< Operator >::BCInterface_Function( const BCInterface_Function& function ) :
        M_parser    ( function.M_parser ),
        M_base      ( function.M_base ),
        M_mapID     ( function.M_mapID )
{
}

// ===================================================
// Methods
// ===================================================
template< typename Operator >
BCInterface_Function< Operator >&
BCInterface_Function< Operator >::operator=( const BCInterface_Function& function )
{
    if ( this != &function )
    {
        M_parser     = function.M_parser;
        M_base       = function.M_base;
        M_mapID      = function.M_mapID;
    }

    return *this;
}

template< typename Operator >
void
BCInterface_Function< Operator >::SetData( const Data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5022 ) << "BCInterface_Function::setData" << "\n";
#endif

    if ( M_parser )
        M_parser->SetString( data.GetBaseString() );
    else
        M_parser.reset( new Parser( data.GetBaseString() ) );

    /*
     * MODE          COMPONENT     FUNCTION      |      COMV.SIZE     ARGUMENTS     INTERFACEFUNCTION
     * ------------------------------------------|---------------------------------------------------
     *                                           |
     * COMPONENT     2             x*y*z         |      1             1             Function
     * FULL          3             x*y*z         |      1             1             Function
     * FULL          1             x*y*z         |      1             1             Function
     * FULL          3             (y*z,x*z,x*y) |      1             3             FunctionID
     * FULL          2             (x,y)         |      1             2             FunctionID
     * COMPONENT     '1 3'         (x,y)         |      2             2             FunctionID
     */

    UInt arguments = M_parser->countSubstring( "," ) + 1;

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface_Function::setFunction            arguments: " << arguments << "\n";
#endif

    if ( arguments == 1 )
        M_base.setFunction( boost::bind( &BCInterface_Function::Function, this, _1, _2, _3, _4, _5 ) );
    else
    {
        //Create the ID map
        if ( data.GetComV().size() > 1 ) // Component
            for ( ID i( 0 ); i < static_cast< ID > ( data.GetComV().size() ); ++i )
                M_mapID[data.GetComV()[i]] = i + 1;
        else
            // if ( data.GetComV().front() == arguments )  Full
            for ( ID i( 1 ); i <= data.GetComV().front(); ++i )
                M_mapID[i] = i;

        M_base.setFunction( boost::bind( &BCInterface_Function::FunctionID, this, _1, _2, _3, _4, _5 ) );
    }
}

// ===================================================
// Get Methods
// ===================================================
template< typename Operator >
typename BCInterface_Function< Operator >::BCFunction_Type&
BCInterface_Function< Operator >::GetBase()
{
    return M_base;
}

template< typename Operator >
Real
BCInterface_Function< Operator >::Function( const Real& t,
                                            const Real& x,
                                            const Real& y,
                                            const Real& z,
                                            const ID& /*id*/)
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface_Function::Function: " << "\n";
    Debug( 5021 ) << "                                                           x: " << x << "\n";
    Debug( 5021 ) << "                                                           y: " << y << "\n";
    Debug( 5021 ) << "                                                           z: " << z << "\n";
    Debug( 5021 ) << "                                                           t: " << t << "\n";
#endif

    M_parser->SetVariable( "t", t );
    M_parser->SetVariable( "x", x );
    M_parser->SetVariable( "y", y );
    M_parser->SetVariable( "z", z );

    this->DataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "                                                evaluate(" << 1 << ") : " << M_parser->Evaluate( 1 ) << "\n";
#endif

    return M_parser->Evaluate( 1 );
}

template< typename Operator >
Real
BCInterface_Function< Operator >::FunctionID( const Real& t,
                                              const Real& x,
                                              const Real& y,
                                              const Real& z,
                                              const ID& id )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface_Function::Function: " << "\n";
    Debug( 5021 ) << "                                                           x: " << x << "\n";
    Debug( 5021 ) << "                                                           y: " << y << "\n";
    Debug( 5021 ) << "                                                           z: " << z << "\n";
    Debug( 5021 ) << "                                                           t: " << t << "\n";
    Debug( 5021 ) << "                                                          id: " << id << "\n";
#endif

    M_parser->SetVariable( "t", t );
    M_parser->SetVariable( "x", x );
    M_parser->SetVariable( "y", y );
    M_parser->SetVariable( "z", z );

    this->DataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "                                                evaluate(" << M_mapID[id] << ") : " << M_parser->Evaluate( M_mapID[id] ) << "\n";
#endif

    return M_parser->Evaluate( M_mapID[id] );
}

} // Namespace LifeV

#endif /* BCInterface_Function_H */
