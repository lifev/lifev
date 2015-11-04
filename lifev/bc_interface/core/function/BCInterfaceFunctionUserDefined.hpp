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
 *  @brief File containing the BCInterfaceFunctionUserDefined class
 *
 *  @date 25-08-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFunctionUserDefined_H
#define BCInterfaceFunctionUserDefined_H 1

// Includes BCInterface classes
#include <lifev/bc_interface/core/function/BCInterfaceFunction.hpp>

namespace LifeV
{

//! BCInterfaceFunctionUserDefined - User defined functions for BCInterface
/*!
 *  @author Cristiano Malossi
 *
 *
 *  The BCInterfaceFunctionUserDefined class provides a set of user defined functions
 *  to be used by the \c BCInterface
 *
 *  <b>DETAILS:</b> <BR>
 *  The constructor of the class takes a string contains the ID of the boundary condition to impose.
 *  The list of available conditions is given through the \c userDefinedFunctions enum. These are:
 *
 *  <ol>
 *      <li> Sin
 *  </ol>
 */
template< typename BcHandlerType, typename PhysicalSolverType >
class BCInterfaceFunctionUserDefined: public virtual BCInterfaceFunction< BcHandlerType, PhysicalSolverType >
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
    typedef BCInterfaceData::parametersContainer_Type                              parametersContainer_Type;

    typedef typename function_Type::bcBase_Type                                    bcBase_Type;

    typedef typename function_Type::data_Type                                      data_Type;
    typedef typename function_Type::dataPtr_Type                                   dataPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    explicit BCInterfaceFunctionUserDefined();

    virtual ~BCInterfaceFunctionUserDefined() {}

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
     * @param timeStep time step
     * @return boundary condition value
     */
    Real functionTime ( const Real& t )
    {
        return functionSelectorTime() ( t );
    };

    //! Function of time and time step
    /*!
     * @param t time
     * @param timeStep time step
     * @return boundary condition value
     */
    Real functionTimeTimeStep ( const Real& t, const Real& timeStep )
    {
        return functionSelectorTimeTimeStep() ( t, timeStep );
    };

    //! Function of time and space
    /*!
     * @param t time
     * @param x x coordinate
     * @param y y coordinate
     * @param z z coordinate
     * @param id id of the boundary condition (not used)
     * @return boundary condition value
     */
    Real functionTimeSpace ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*id*/)
    {
        return functionSelectorTimeSpaceID() ( t, x, y, z, 0 );
    };

    //! Function of time and space with ID
    /*!
     * @param t time
     * @param x x coordinate
     * @param y y coordinate
     * @param z z coordinate
     * @param id id of the boundary condition
     * @return boundary condition value
     */
    Real functionTimeSpaceID ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id )
    {
        return functionSelectorTimeSpaceID() ( t, x, y, z, id );
    };


    //@}


    //! @name Set Methods
    //@{

    //! Set data for boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setData ( const dataPtr_Type& data );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionUserDefined ( const BCInterfaceFunctionUserDefined& function);

    BCInterfaceFunctionUserDefined& operator= ( const BCInterfaceFunctionUserDefined& function );

    //@}


    //! @name Private Methods
    //@{

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

    //! Get the selected function of time and space with ID
    /*!
     * @return boundary function
     */
    boundaryFunctionTimeSpaceID_Type functionSelectorTimeSpaceID();

    //! Sinusoidal function
    /*!
     * @param t time
     * @param x x coordinate
     * @param y y coordinate
     * @param z z coordinate
     * @param id id of the boundary condition (not used)
     * @return boundary condition value
     */
    Real functionSin ( const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/ );

    //@}

    enum userDefinedFunctions
    {
        Sin
    };

    userDefinedFunctions     M_functionType;

    parametersContainer_Type M_parameters;

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename BcHandlerType, typename PhysicalSolverType >
inline BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >* createBCInterfaceFunctionUserDefined()
{
    return new BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType > ();
}

// ===================================================
// Constructor
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::BCInterfaceFunctionUserDefined() :
    function_Type   (),
    M_functionType  (),
    M_parameters    ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5026 ) << "BCInterfaceFunctionUserDefined::BCInterfaceFunctionUserDefined()" << "\n";
#endif

}

// ===================================================
// Private methods
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
void
BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::setData ( const dataPtr_Type& data )
{
    //Set mapFunction
    std::map< std::string, userDefinedFunctions > mapFunction;
    mapFunction["Sin"] = Sin;

    // Retrieving the strings
    M_functionType = mapFunction[ data->baseString() ];

    // Set parameters
    M_parameters = data->parameters();
}

template< typename BcHandlerType, typename PhysicalSolverType >
typename BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::boundaryFunctionTime_Type
BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::functionSelectorTime()
{
    switch ( M_functionType )
    {
        case Sin:

            return std::bind ( &BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::functionSin, this, std::placeholders::_1, std::placeholders::_1, std::placeholders::_1, std::placeholders::_1, std::placeholders::_1 );

        default:

            std::cout << " !!! Warning: user defined function " << M_functionType << " not assigned !!!" << std::endl;

            return 0;
    }
}

template< typename BcHandlerType, typename PhysicalSolverType >
typename BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::boundaryFunctionTimeTimeStep_Type
BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::functionSelectorTimeTimeStep()
{
    switch ( M_functionType )
    {
        case Sin:

            return std::bind ( &BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::functionSin, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_1, std::placeholders::_1, std::placeholders::_1 );

        default:

            std::cout << " !!! Warning: user defined function " << M_functionType << " not assigned !!!" << std::endl;

            return 0;
    }
}

template< typename BcHandlerType, typename PhysicalSolverType >
typename BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::boundaryFunctionTimeSpaceID_Type
BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::functionSelectorTimeSpaceID()
{
    switch ( M_functionType )
    {
        case Sin:

            return std::bind ( &BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::functionSin, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );

        default:

            std::cout << " !!! Warning: user defined function " << M_functionType << " not assigned !!!" << std::endl;

            return 0;
    }
}

template< typename BcHandlerType, typename PhysicalSolverType >
Real
BCInterfaceFunctionUserDefined< BcHandlerType, PhysicalSolverType >::functionSin ( const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return M_parameters[0] + M_parameters[1] * std::sin ( M_parameters[2] + 2 * M_PI * t / M_parameters[3] );
}

} // Namespace LifeV

#endif /* BCInterfaceFunctionUserDefined_H */
