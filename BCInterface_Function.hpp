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
 *  @brief File containing the BCInterface_Function class
 *
 *  @date 06-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
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
class BCInterface_Function
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterface_Data                                                      data_Type;
    typedef BCFunctionBase                                                        bcFunction_Type;
    typedef parser::Parser                                                        parser_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    BCInterface_Function();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface_Function( const data_Type& data );

    //! Destructor
    virtual ~BCInterface_Function() {}

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

    BCInterface_Function( const BCInterface_Function& function );

    BCInterface_Function& operator=( const BCInterface_Function& function );

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
inline BCInterface_Function< PhysicalSolverType >* createBCInterface_Function()
{
    return new BCInterface_Function< PhysicalSolverType > ();
}

// ===================================================
// Constructor
// ===================================================
template< typename PhysicalSolverType >
BCInterface_Function< PhysicalSolverType >::BCInterface_Function() :
        M_parser    (),
        M_base      (),
        M_mapID     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface_Function::BCInterface_Function()" << "\n";
#endif

}

template< typename PhysicalSolverType >
BCInterface_Function< PhysicalSolverType >::BCInterface_Function( const data_Type& data ) :
        M_parser    (),
        M_base      (),
        M_mapID     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface_Function::BCInterface_Function( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Set Methods
// ===================================================
template< typename PhysicalSolverType >
void
BCInterface_Function< PhysicalSolverType >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5022 ) << "BCInterface_Function::setData" << "\n";
#endif

    if ( M_parser )
        M_parser->SetString( data.baseString() );
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
    Debug( 5021 ) << "BCInterface_Function::setFunction            arguments: " << arguments << "\n";
#endif

    if ( arguments == 1 )
        M_base.setFunction( boost::bind( &BCInterface_Function::function, this, _1, _2, _3, _4, _5 ) );
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

        M_base.setFunction( boost::bind( &BCInterface_Function::functionID, this, _1, _2, _3, _4, _5 ) );
    }
}

// ===================================================
// Get Methods
// ===================================================
template< typename PhysicalSolverType >
Real
BCInterface_Function< PhysicalSolverType >::function( const Real& t,
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

    this->dataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "                                                evaluate(" << 1 << ") : " << M_parser->Evaluate( 1 ) << "\n";
#endif

    return M_parser->Evaluate( 1 );
}

template< typename PhysicalSolverType >
Real
BCInterface_Function< PhysicalSolverType >::functionID( const Real& t,
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

    this->dataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "                                                evaluate(" << M_mapID[id] << ") : " << M_parser->Evaluate( M_mapID[id] ) << "\n";
#endif

    return M_parser->Evaluate( M_mapID[id] );
}

} // Namespace LifeV

#endif /* BCInterface_Function_H */
