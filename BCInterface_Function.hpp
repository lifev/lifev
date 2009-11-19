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

#include <lifemc/lifefem/BCInterface_Definitions.hpp>
#include <lifemc/lifefem/BCInterface_Data.hpp>

#include <lifemc/lifecore/SpiritParser.hpp>

namespace LifeV {

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
template< typename Operator >
class BCInterface_Function
{
public:

    typedef boost::function< Real( Real const& t,
                                   Real const& x,
                                   Real const& y,
                                   Real const& z,
                                   ID const& id ) > function_type;

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    BCInterface_Function();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface_Function( const BCInterface_Data< Operator >& data );

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
    virtual void SetData( const BCInterface_Data< Operator >& data );

    //! Compare function
    /*!
     * @param data BC data loaded from GetPot file
     * @return true if the functions are equal, false if they aren't
     */
    virtual bool Compare( const BCInterface_Data< Operator >& data );

    //@}


    //! @name Get functions
    //@{

    //! Get the base of the boundary condition
    BCFunctionBase& GetBase()
    {
        return M_base;
    }

    //@}

protected:

    std::string                       M_baseString;
    BCComV                            M_comV;
    boost::shared_ptr< SpiritParser > M_parser;

    //! @name Protected functions
    //@{

    //! dataInterpolation
    virtual inline void DataInterpolation() {}

    //@}

private:

    BCFunctionBase                   M_base;
    std::map< ID, ID >               M_mapID;

    //! @name Private functions
    //@{

    //! SetFunction
    inline void SetFunction();

    //! Function
    Real Function( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*id*/);

    //! FunctionID
    Real FunctionID( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id );

    //@}
};

//! Factory create function
template< typename Operator >
inline BCInterface_Function< Operator >* createFunction()
{
    return new BCInterface_Function< Operator > ();
}

// ===================================================
// Constructor
// ===================================================
template< typename Operator >
BCInterface_Function< Operator >::BCInterface_Function() :
    M_baseString(),
    M_comV      (),
    M_parser    (),
    M_base      (),
    M_mapID     ()
{

#ifdef DEBUG
    Debug( 5021 ) << "BCInterface_Function::BCInterface_Function()" << "\n";
#endif

}

template< typename Operator >
BCInterface_Function< Operator >::BCInterface_Function( const BCInterface_Data< Operator >& data ) :
    M_baseString(),
    M_comV      (),
    M_parser    (),
    M_base      (),
    M_mapID     ()
{

#ifdef DEBUG
    Debug( 5021 ) << "BCInterface_Function::BCInterface_Function( data )" << "\n";
#endif

    this->SetData( data );
}

template< typename Operator >
BCInterface_Function< Operator >::BCInterface_Function( const BCInterface_Function& function ) :
    M_baseString( function.M_baseString ),
    M_comV      ( function.M_comV ),
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
        M_baseString = function.M_baseString;
        M_comV       = function.M_comV;
        M_parser     = function.M_parser;
        M_base       = function.M_base;
        M_mapID      = function.M_mapID;
    }

    return *this;
}

template< typename Operator >
void
BCInterface_Function< Operator >::SetData( const BCInterface_Data< Operator >& data )
{

#ifdef DEBUG
    Debug( 5022 ) << "BCInterface_Function::setData" << "\n";
#endif

    M_comV = data.GetComV();
    M_baseString = data.GetBaseString();

    //boost::shared_ptr<SpiritParser> emptyParser( );
    //if ( M_parser == emptyParser )
    if ( M_parser )
        M_parser->setString( M_baseString );
    else
        M_parser.reset( new SpiritParser( M_baseString ) );

    SetFunction();
}

template< typename Operator >
bool
BCInterface_Function< Operator >::Compare( const BCInterface_Data< Operator >& data )
{
    return M_baseString.compare( data.GetBaseString() ) == 0 && M_comV == data.GetComV();
}

// ===================================================
// Private functions
// ===================================================
template< typename Operator >
inline void BCInterface_Function< Operator >::SetFunction()
{
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

#ifdef DEBUG
    Debug( 5021 ) << "BCInterface_Function::setFunction            arguments: " << arguments << "\n";
#endif

    if ( arguments == 1 )
        M_base.setFunction( boost::bind( &BCInterface_Function::Function, this, _1, _2, _3, _4, _5 ) );
    else
    {
        //Create the ID map
        if ( M_comV.size() > 1 ) // Component
            for ( ID i( 0 ); i < static_cast< ID > ( M_comV.size() ); ++i )
                M_mapID[M_comV[i]] = i + 1;
        else
            // if ( M_comV.front() == arguments )  Full
            for ( ID i( 1 ); i <= M_comV.front(); ++i )
                M_mapID[i] = i;

        M_base.setFunction( boost::bind( &BCInterface_Function::FunctionID, this, _1, _2, _3, _4, _5 ) );
    }
}

template< typename Operator >
Real
BCInterface_Function< Operator >::Function( const Real& t,
                                                const Real& x,
                                                const Real& y,
                                                const Real& z,
                                                const ID& /*id*/)
{

#ifdef DEBUG
    Debug( 5021 ) << "BCInterface_Function::Function: " << "\n";
    Debug( 5021 ) << "                                                           x: " << x << "\n";
    Debug( 5021 ) << "                                                           y: " << y << "\n";
    Debug( 5021 ) << "                                                           z: " << z << "\n";
    Debug( 5021 ) << "                                                           t: " << t << "\n";
#endif

    M_parser->setVariable( "t", t );
    M_parser->setVariable( "x", x );
    M_parser->setVariable( "y", y );
    M_parser->setVariable( "z", z );

    this->DataInterpolation();

#ifdef DEBUG
    Debug( 5021 ) << "                                                evaluate(" << 1 << ") : " << M_parser->evaluate( 1 ) << "\n";
#endif

    return M_parser->evaluate( 1 );
}

template< typename Operator >
Real
BCInterface_Function< Operator >::FunctionID( const Real& t,
                                                  const Real& x,
                                                  const Real& y,
                                                  const Real& z,
                                                  const ID& id )
{

#ifdef DEBUG
    Debug( 5021 ) << "BCInterface_Function::Function: " << "\n";
    Debug( 5021 ) << "                                                           x: " << x << "\n";
    Debug( 5021 ) << "                                                           y: " << y << "\n";
    Debug( 5021 ) << "                                                           z: " << z << "\n";
    Debug( 5021 ) << "                                                           t: " << t << "\n";
    Debug( 5021 ) << "                                                          id: " << id << "\n";
#endif

    M_parser->setVariable( "t", t );
    M_parser->setVariable( "x", x );
    M_parser->setVariable( "y", y );
    M_parser->setVariable( "z", z );

    this->DataInterpolation();

#ifdef DEBUG
    Debug( 5021 ) << "                                                evaluate(" << M_mapID[id] << ") : " << M_parser->evaluate( M_mapID[id] ) << "\n";
#endif

    return M_parser->evaluate( M_mapID[id] );
}

} // Namespace LifeV

#endif /* BCInterface_Function_H */
