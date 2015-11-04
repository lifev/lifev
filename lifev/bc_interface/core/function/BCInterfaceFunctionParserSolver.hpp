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
 *  @brief File containing the BCInterfaceFunctionParserSolver class
 *
 *  @date 24-08-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFunctionParserSolver_H
#define BCInterfaceFunctionParserSolver_H 1

// BCInterface includes
#include <lifev/bc_interface/core/function/BCInterfaceFunctionParser.hpp>

namespace LifeV
{

//! BCInterfaceFunctionParserSolver - LifeV boundary condition function file wrapper for \c BCInterface
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between the \c BCInterface, the \c Parser, and a general
 *  LifeV physical solver (such as \c OseenSolver or \c FSISolver). It allows to construct LifeV
 *  function types for boundary conditions, using a functions string loaded from
 *  a \c GetPot file in which are present some physical solver parameters.
 *
 *  The class can be used in two ways:
 *
 *  <ol>
 *      <li> first hereditate it and then implement the template specialization for the methods \c createAccessList() and \c updatePhysicalSolverVariables();
 *      <li> manually setting the variables by using the \c setVariable() method.
 *  </ol>
 *
 *  See \c BCInterfaceFunctionParser class for more details.
 *
 *  <b>AVAILABLE VARIABLES</b> <BR>
 *  Current available variables are:
 *
 *  <ul>
 *      <li> f_timeStep
 *        <li> f_area
 *        <li> f_density
 *        <li> f_flux
 *        <li> f_pressure
 *        <li> f_viscosity
 *        <li> f_venousPressure
 *        <li> s_density
 *        <li> s_poisson
 *        <li> s_thickness
 *        <li> s_young
 *        <li> s_externalPressure
 *    </ul>
 *
 *    Of course, some of those variables are available only for fluid problems, other only for solid problems.
 */
template< typename BcHandlerType, typename PhysicalSolverType >
class BCInterfaceFunctionParserSolver: public virtual BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef BcHandlerType                                                          bcHandler_Type;
    typedef PhysicalSolverType                                                     physicalSolver_Type;

    typedef std::shared_ptr< physicalSolver_Type >                               physicalSolverPtr_Type;
    typedef BCInterfaceFunction< bcHandler_Type, physicalSolver_Type >             function_Type;
    typedef BCInterfaceFunctionParser< bcHandler_Type, physicalSolver_Type >       functionParser_Type;
    typedef typename PhysicalSolverType::solutionPtr_Type                          solutionPtr_Type;

    typedef typename function_Type::data_Type                                      data_Type;
    typedef typename function_Type::dataPtr_Type                                   dataPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceFunctionParserSolver();

    //! Destructor
    virtual ~BCInterfaceFunctionParserSolver() {}

    //@}


    //! @name Methods
    //@{

    //! Update the solver variables
    /*!
     *  <b>NOTE:</b> A template specialization of this method should be provided for each solver.
     */
    void updatePhysicalSolverVariables()
    {
        std::cout << " !!! WARNING: updatePhysicalSolverVariables() is not defined for the selected solver. !!!" << std::endl;
    }

    //@}


    //! @name Set Methods
    //@{

    //! Set data for boundary conditions
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void setData ( const dataPtr_Type& data );

    //! Set the physical solver
    /*!
     * @param physicalSolver physical solver
     */
    void setPhysicalSolver ( const physicalSolverPtr_Type& physicalSolver )
    {
        M_physicalSolver = physicalSolver;
    }

    //! Set solution
    /*!
     * @param solution The solution container of the 1D problem
     */
    void setSolution ( const solutionPtr_Type& solution )
    {
        M_solution = solution;
    }

    //! Set variable function
    /*!
     * @param name name of the variable
     * @param value value of the variable
     */
    void setVariable ( const std::string& name, const Real& value )
    {
        functionParser_Type::M_parser->setVariable ( name, value );
    }

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Create the list of variables for the physical solver.
    /*!
     *  NOTE: A template specialization of this method should be provided for each solver.
     */
    void createAccessList ( const dataPtr_Type& /*data*/ )
    {
        std::cout << " !!! WARNING: createAccessList() is not defined for the selected solver. !!!" << std::endl;
    }

    //@}

    //List of all available operators (f: fluid; s: solid;)
    enum physicalSolverList
    {
        f_timeStep,
        f_area,
        f_flux,
        f_density,
        f_pressure,
        f_viscosity,
        f_venousPressure,
        s_density,
        s_poisson,
        s_thickness,
        s_young,
        s_externalPressure
    };

    physicalSolverPtr_Type                     M_physicalSolver;
    solutionPtr_Type                           M_solution;

    ID                                         M_boundaryID;
    std::set< physicalSolverList >             M_list;

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionParserSolver ( const BCInterfaceFunctionParserSolver& function );

    BCInterfaceFunctionParserSolver& operator= ( const BCInterfaceFunctionParserSolver& function );

    //@}


    //! @name Private Methods
    //@{

    void createFluidMap ( std::map< std::string, physicalSolverList >& mapList );
    void createSolidMap ( std::map< std::string, physicalSolverList >& mapList );
    void createList ( const std::map< std::string, physicalSolverList >& mapList, const dataPtr_Type& data );

    void switchErrorMessage ( const std::string& operatorType )
    {
        std::cout << "ERROR: Invalid variable type for " << operatorType << " FunctionSolver" << std::endl;
    }

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename BcHandlerType, typename PhysicalSolverType >
inline BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >* createBCInterfaceFunctionParserSolver()
{
    return new BCInterfaceFunctionParserSolver< BcHandlerType, PhysicalSolverType > ();
}

// ===================================================
// Constructors
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
BCInterfaceFunctionParserSolver< BcHandlerType, PhysicalSolverType >::BCInterfaceFunctionParserSolver() :
    function_Type                    (),
    functionParser_Type              (),
    M_physicalSolver                 (),
    M_solution                       (),
    M_boundaryID                     (),
    M_list                           ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver::BCInterfaceFunctionSolver()" << "\n";
#endif

}



// ===================================================
// Set Methods
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
void
BCInterfaceFunctionParserSolver< BcHandlerType, PhysicalSolverType >::setData ( const dataPtr_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver::setData( data )" << "\n";
#endif

    M_boundaryID = data->boundaryID();

    functionParser_Type::setData ( data );

    createAccessList ( data );
}



// ===================================================
// Private Methods
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
inline void
BCInterfaceFunctionParserSolver< BcHandlerType, PhysicalSolverType >::createFluidMap ( std::map< std::string, physicalSolverList >& mapList )
{
    mapList["f_timeStep"]       = f_timeStep;
    mapList["f_area"]           = f_area;
    mapList["f_density"]        = f_density;
    mapList["f_flux"]           = f_flux;
    mapList["f_pressure"]       = f_pressure;
    mapList["f_viscosity"]      = f_viscosity;
    mapList["f_venousPressure"] = f_venousPressure;
}

template< typename BcHandlerType, typename PhysicalSolverType >
inline void
BCInterfaceFunctionParserSolver< BcHandlerType, PhysicalSolverType >::createSolidMap ( std::map< std::string, physicalSolverList >& mapList )
{
    mapList["s_density"]          = s_density;
    mapList["s_poisson"]          = s_poisson;
    mapList["s_thickness"]        = s_thickness;
    mapList["s_young"]            = s_young;
    mapList["s_externalPressure"] = s_externalPressure;
}

template< typename BcHandlerType, typename PhysicalSolverType >
inline void
BCInterfaceFunctionParserSolver< BcHandlerType, PhysicalSolverType >::createList ( const std::map< std::string, physicalSolverList >& mapList, const dataPtr_Type& data )
{
    M_list.clear();
    for ( typename std::map< std::string, physicalSolverList >::const_iterator j = mapList.begin(); j != mapList.end(); ++j )
        if ( boost::find_first ( data->baseString(), j->first ) )
        {
            M_list.insert ( j->second );
        }
}

} // Namespace LifeV

#endif /* BCInterfaceFunctionParserSolver_H */
