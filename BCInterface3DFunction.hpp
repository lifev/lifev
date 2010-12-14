//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing the BCInterface3DFunction class
 *
 *  @date 06-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface3DFunction_H
#define BCInterface3DFunction_H 1

#include <lifemc/lifesolver/BCInterface3DDefinitions.hpp>
#include <lifemc/lifesolver/BCInterface3DData.hpp>

#include <lifemc/lifecore/Parser.hpp>

namespace LifeV
{

//! BCInterface3DFunction - LifeV bcFunction wrapper for BCInterface
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface and grammar parser. It allows to construct LifeV
 *  functions type for boundary conditions, using a functions string loaded from a GetPot file.
 *
 *  <b>DETAILS:</b>
 *
 *  By default the stringSeparator is set to semicolon ";".
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
template< typename PhysicalSolverType >
class BCInterface3DFunction
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterface3DData                                                     data_Type;
    typedef BCFunctionBase                                                        bcFunction_Type;
    typedef Parser                                                                parser_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    explicit BCInterface3DFunction();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    explicit BCInterface3DFunction( const data_Type& data );

    //! Destructor
    virtual ~BCInterface3DFunction() {}

    //@}


    //! @name Set Methods
    //@{

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void setData( const data_Type& data );

    //@}


    //! @name Get Methods
    //@{

    //! Get the base of the boundary condition
    /*!
     * @return boundary condition base
     */
    bcFunction_Type& base() { return M_base; }

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! dataInterpolation
    virtual void dataInterpolation() {}

    //@}

    boost::shared_ptr< parser_Type > M_parser;

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface3DFunction( const BCInterface3DFunction& function );

    BCInterface3DFunction& operator=( const BCInterface3DFunction& function );

    //@}


    //! @name Private Methods
    //@{

    //! function
    Real function( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*id*/);

    //! functionID
    Real functionID( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id );

    //@}

    bcFunction_Type                  M_base;
    std::map< ID, ID >               M_mapID;

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename PhysicalSolverType >
inline BCInterface3DFunction< PhysicalSolverType >* createBCInterface3DFunction()
{
    return new BCInterface3DFunction< PhysicalSolverType > ();
}

// ===================================================
// Constructor
// ===================================================
template< typename PhysicalSolverType >
BCInterface3DFunction< PhysicalSolverType >::BCInterface3DFunction() :
        M_parser    (),
        M_base      (),
        M_mapID     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface3DFunction::BCInterface3DFunction()" << "\n";
#endif

}

template< typename PhysicalSolverType >
BCInterface3DFunction< PhysicalSolverType >::BCInterface3DFunction( const data_Type& data ) :
        M_parser    (),
        M_base      (),
        M_mapID     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface3DFunction::BCInterface3DFunction( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Set Methods
// ===================================================
template< typename PhysicalSolverType >
void
BCInterface3DFunction< PhysicalSolverType >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5022 ) << "BCInterface3DFunction::setData" << "\n";
#endif

    if ( M_parser )
        M_parser->setString( data.baseString() );
    else
        M_parser.reset( new parser_Type( data.baseString() ) );

    /*
     * MODE          COMPONENT     FUNCTION      |      COMV.SIZE     ARGUMENTS     INTERFACEFUNCTION
     * ------------------------------------------|---------------------------------------------------
     *                                           |
     * COMPONENT     2             x*y*z         |      1             1             function
     * FULL          3             x*y*z         |      1             1             function
     * FULL          1             x*y*z         |      1             1             function
     * FULL          3             (y*z,x*z,x*y) |      1             3             functionID
     * FULL          2             (x,y)         |      1             2             functionID
     * COMPONENT     '1 3'         (x,y)         |      2             2             functionID
     */

    UInt arguments = M_parser->countSubstring( "," ) + 1;

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface3DFunction::setFunction            arguments: " << arguments << "\n";
#endif

    if ( arguments == 1 )
        M_base.setFunction( boost::bind( &BCInterface3DFunction::function, this, _1, _2, _3, _4, _5 ) );
    else
    {
        //Create the ID map
        if ( data.comV().size() > 1 ) // Component
            for ( ID i( 0 ); i < static_cast< ID > ( data.comV().size() ); ++i )
                M_mapID[data.comV()[i]] = i + 1;
        else
            // if ( data.comV().front() == arguments )  Full
            for ( ID i( 1 ); i <= data.comV().front(); ++i )
                M_mapID[i] = i;

        M_base.setFunction( boost::bind( &BCInterface3DFunction::functionID, this, _1, _2, _3, _4, _5 ) );
    }
}

// ===================================================
// Private Methods
// ===================================================
template< typename PhysicalSolverType >
Real
BCInterface3DFunction< PhysicalSolverType >::function( const Real& t,
                                            const Real& x,
                                            const Real& y,
                                            const Real& z,
                                            const ID& /*id*/)
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface3DFunction::Function: " << "\n";
    Debug( 5021 ) << "                                                           x: " << x << "\n";
    Debug( 5021 ) << "                                                           y: " << y << "\n";
    Debug( 5021 ) << "                                                           z: " << z << "\n";
    Debug( 5021 ) << "                                                           t: " << t << "\n";
#endif

    M_parser->setVariable( "t", t );
    M_parser->setVariable( "x", x );
    M_parser->setVariable( "y", y );
    M_parser->setVariable( "z", z );

    this->dataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "                                                evaluate(" << 1 << ") : " << M_parser->evaluate( 1 ) << "\n";
#endif

    return M_parser->evaluate( 1 );
}

template< typename PhysicalSolverType >
Real
BCInterface3DFunction< PhysicalSolverType >::functionID( const Real& t,
                                              const Real& x,
                                              const Real& y,
                                              const Real& z,
                                              const ID& id )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface3DFunction::Function: " << "\n";
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

    this->dataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "                                                evaluate(" << M_mapID[id] << ") : " << M_parser->evaluate( M_mapID[id] ) << "\n";
#endif

    return M_parser->evaluate( M_mapID[id] );
}

} // Namespace LifeV

#endif /* BCInterface3DFunction_H */
