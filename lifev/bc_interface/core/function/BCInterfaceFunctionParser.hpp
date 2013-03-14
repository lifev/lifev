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
 *  @brief File containing the BCInterfaceFunctionParser class
 *
 *  @date 06-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFunctionParser_H
#define BCInterfaceFunctionParser_H 1

// Includes parser classes
#include <lifev/core/util/Parser.hpp>

// Includes BCInterface classes
#include <lifev/bc_interface/core/function/BCInterfaceFunction.hpp>

namespace LifeV
{

//! BCInterfaceFunctionParser - LifeV boundary condition function wrapper for \c BCInterface
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between the \c BCInterface and the \c Parser. It allows to construct LifeV
 *  function types for boundary conditions, using a function string loaded from a \c GetPot file.
 *
 *  <b>DETAILS</b> <BR>
 *  By default the string separator is set to semicolon ";".
 *
 *  The function string has to be in this form:
 *
 *  <CODE>
 *  function = '[u, v, w]'
 *  </CODE>
 *
 *  where \f$ u=u(x,y,z,t) \f$, \f$ v=v(x,y,z,t) \f$,  and \f$ w=w(x,y,z,t) \f$. They are separated by commas,
 *  as shown in the following example
 *
 *  <CODE>
 *  function = '[x^2 + y^2, 0, 2*sin(2*pi*t)]'
 *  </CODE>
 *
 *  It is possible to define constants separated by the string separator before the expression
 *
 *  <CODE>
 *  function = 'a=5.67436; [x^2+y^2,0,a*sin(2*pi*t)^a]'
 *  </CODE>
 *
 *  <b>NOTE</b> <BR>
 *  In the boundary condition file, if you have three component with the same expression
 *  (the same function) you can both write:
 *
 *  <CODE>
 *  function = '[0, 0, 0]'
 *  </CODE>
 *
 *  or
 *
 *  <CODE>
 *  function = 0
 *  </CODE>
 *
 *  However the second way is more efficient during execution.
 */
template< typename BcHandlerType, typename PhysicalSolverType >
class BCInterfaceFunctionParser: public virtual BCInterfaceFunction< BcHandlerType, PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef BcHandlerType                                                          bcHandler_Type;
    typedef PhysicalSolverType                                                     physicalSolver_Type;

    typedef BCInterfaceFunction< bcHandler_Type, physicalSolver_Type >             function_Type;
    typedef typename function_Type::boundaryFunctionTime_Type                      boundaryFunctionTime_Type;
    typedef typename function_Type::boundaryFunctionTimeTimeStep_Type              boundaryFunctionTimeTimeStep_Type;
    typedef typename function_Type::boundaryFunctionTimeSpaceID_Type               boundaryFunctionTimeSpaceID_Type;
    typedef Parser                                                                 parser_Type;
    typedef boost::shared_ptr< parser_Type >                                       parserPtr_Type;

    typedef typename function_Type::data_Type                                      data_Type;
    typedef typename function_Type::dataPtr_Type                                   dataPtr_Type;

    typedef typename function_Type::bcBase_Type                                    bcBase_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    explicit BCInterfaceFunctionParser();

    //! Destructor
    virtual ~BCInterfaceFunctionParser() {}

    //@}


    //! @name Methods
    //@{

    //! Assign the function to the base of the \c BCHandler
    /*!
     * @param base base of the boundary condition
     */
    void assignFunction ( bcBase_Type& base );

    //! Function of time
    /*!
     * @param t time
     * @return boundary condition value
     */
    Real functionTime ( const Real& t );

    //! Function of time and time step
    /*!
     * @param t time
     * @param timeStep time step
     * @return boundary condition value
     */
    Real functionTimeTimeStep ( const Real& t, const Real& timeStep );

    //! Function of time and space
    /*!
     * @param t time
     * @param x x coordinate
     * @param y y coordinate
     * @param z z coordinate
     * @param id id of the boundary condition (not used)
     * @return boundary condition value
     */
    Real functionTimeSpace ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*id*/);

    //! Function of time and space with ID
    /*!
     * @param t time
     * @param x x coordinate
     * @param y y coordinate
     * @param z z coordinate
     * @param id id of the boundary condition
     * @return boundary condition value
     */
    Real functionTimeSpaceID ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id );

    //@}


    //! @name Set Methods
    //@{

    //! Set data for boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    virtual void setData ( const dataPtr_Type& data );

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! dataInterpolation
    virtual void dataInterpolation() {}

    //@}

    parserPtr_Type M_parser;

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionParser ( const BCInterfaceFunctionParser& function );

    BCInterfaceFunctionParser& operator= ( const BCInterfaceFunctionParser& function );

    //@}

    //! @name Private Methods
    //@{

    //! Setup the parser using the data container
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setupParser ( const dataPtr_Type& data );

    //! Get the selected function of time
    /*!
     * @return boundary function
     */
    boundaryFunctionTime_Type functionSelectorTime();

    //! Get the selected function of time and time step
    /*!
     * @return boundary function
     */
    boundaryFunctionTimeTimeStep_Type functionSelectorTimeTimeStep();

    //! Get the selected function of time space and ID.
    /*!
     * @return boundary function
     */
    boundaryFunctionTimeSpaceID_Type functionSelectorTimeSpaceID();

    //@}

    std::map< ID, ID >               M_mapID;

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename BcHandlerType, typename PhysicalSolverType >
inline BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >* createBCInterfaceFunctionParser()
{
    return new BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType > ();
}

// ===================================================
// Constructor
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::BCInterfaceFunctionParser() :
    function_Type   (),
    M_parser        (),
    M_mapID         ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5021 ) << "BCInterfaceFunction::BCInterfaceFunction()" << "\n";
#endif

}

// ===================================================
// Methods
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
Real
BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::functionTime ( const Real& t )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5021 ) << "BCInterfaceFunction::functionTime: " << "\n";
    debugStream ( 5021 ) << "                                                           t: " << t << "\n";
#endif

    M_parser->setVariable ( "t", t );

    this->dataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5021 ) << "                                                evaluate( 0 ) : " << M_parser->evaluate ( 0 ) << "\n";
#endif

    return M_parser->evaluate ( 0 );
}

template< typename BcHandlerType, typename PhysicalSolverType >
Real
BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::functionTimeTimeStep ( const Real& t, const Real& timeStep )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5021 ) << "BCInterfaceFunction::functionTime: " << "\n";
    debugStream ( 5021 ) << "                                                           t: " << t << "\n";
    debugStream ( 5021 ) << "                                                           timeStep: " << timeStep << "\n";
#endif

    M_parser->setVariable ( "t", t );
    M_parser->setVariable ( "timeStep", timeStep );

    this->dataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5021 ) << "                                                evaluate( 0 ) : " << M_parser->evaluate ( 0 ) << "\n";
#endif

    return M_parser->evaluate ( 0 );
}

template< typename BcHandlerType, typename PhysicalSolverType >
Real
BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::functionTimeSpace ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*id*/)
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5021 ) << "BCInterfaceFunction::functionTimeSpace: " << "\n";
    debugStream ( 5021 ) << "                                                           x: " << x << "\n";
    debugStream ( 5021 ) << "                                                           y: " << y << "\n";
    debugStream ( 5021 ) << "                                                           z: " << z << "\n";
    debugStream ( 5021 ) << "                                                           t: " << t << "\n";
#endif

    M_parser->setVariable ( "t", t );
    M_parser->setVariable ( "x", x );
    M_parser->setVariable ( "y", y );
    M_parser->setVariable ( "z", z );

    this->dataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5021 ) << "                                                evaluate( 0 ) : " << M_parser->evaluate ( 0 ) << "\n";
#endif

    return M_parser->evaluate ( 0 );
}

template< typename BcHandlerType, typename PhysicalSolverType >
Real
BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::functionTimeSpaceID ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5021 ) << "BCInterfaceFunction::functionTimeSpaceID: " << "\n";
    debugStream ( 5021 ) << "                                                           x: " << x << "\n";
    debugStream ( 5021 ) << "                                                           y: " << y << "\n";
    debugStream ( 5021 ) << "                                                           z: " << z << "\n";
    debugStream ( 5021 ) << "                                                           t: " << t << "\n";
    debugStream ( 5021 ) << "                                                          id: " << id << "\n";
#endif

    M_parser->setVariable ( "t", t );
    M_parser->setVariable ( "x", x );
    M_parser->setVariable ( "y", y );
    M_parser->setVariable ( "z", z );

    this->dataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5021 ) << "                                                evaluate(" << M_mapID[id] << ") : " << M_parser->evaluate ( M_mapID[id] ) << "\n";
#endif

    return M_parser->evaluate ( M_mapID[id] );
}

// ===================================================
// Private Methods
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
void
BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::setupParser ( const dataPtr_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5022 ) << "BCInterfaceFunction::setParser" << "\n";
#endif

    if ( M_parser )
    {
        M_parser->setString ( data->baseString() );
    }
    else
    {
        M_parser.reset ( new parser_Type ( data->baseString() ) );
    }
}

template< typename BcHandlerType, typename PhysicalSolverType >
typename BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::boundaryFunctionTime_Type
BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::functionSelectorTime()
{
    return boost::bind ( &BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::functionTime, this, _1 );
}

template< typename BcHandlerType, typename PhysicalSolverType >
typename BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::boundaryFunctionTimeTimeStep_Type
BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::functionSelectorTimeTimeStep()
{
    return boost::bind ( &BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::functionTimeTimeStep, this, _1, _2 );
}

template< typename BcHandlerType, typename PhysicalSolverType >
typename BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::boundaryFunctionTimeSpaceID_Type
BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::functionSelectorTimeSpaceID()
{
    if ( M_parser->countSubstring ( "," ) )
    {
        return boost::bind ( &BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::functionTimeSpaceID, this, _1, _2, _3, _4, _5 );
    }
    else
    {
        return boost::bind ( &BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >::functionTimeSpace, this, _1, _2, _3, _4, _5 );
    }
}

} // Namespace LifeV

#endif /* BCInterfaceFunctionParser_H */
