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

// Oseen includes
#include <lifev/navier_stokes/solver/OseenSolverShapeDerivative.hpp>

// FSI includes
#include <lifev/fsi/solver/FSIOperator.hpp>

// OneDFSI includes
#include <lifev/one_d_fsi/solver/OneDFSISolver.hpp>

// BCInterface includes
#include <lifev/bc_interface/function/BCInterfaceFunctionParser.hpp>

namespace LifeV
{

#ifndef MULTISCALE_IS_IN_LIFEV
//! ZeroDimensionalData - Temporary data container for ZeroDimensionalModels until they are defined only in LifeV
/*!
 * This class will be removed when the Multiscale framework will be ported in LifeV
 */
class ZeroDimensionalTemporaryData
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalTemporaryData() : M_fluidVenousPressure() {}

    //! Destructor
    virtual ~ZeroDimensionalTemporaryData() {}

    //@}


    //! @name Set Methods
    //@{

    //! Set the global fluid venous pressure.
    /*!
     * @return venous pressure of the fluid.
     */
    void setFluidVenousPressure( const Real& fluidVenousPressure ) { M_fluidVenousPressure = fluidVenousPressure; }

    //@}


    //! @name Get Methods
    //@{

    //! Get the global fluid venous pressure.
    /*!
     * @return venous pressure of the fluid.
     */
    const Real& fluidVenousPressure() const { return M_fluidVenousPressure; }

    //@}

private:

    Real                                M_fluidVenousPressure;
};
#endif

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
template< class PhysicalSolverType >
class BCInterfaceFunctionParserSolver: public virtual BCInterfaceFunctionParser< PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterfaceFunction< physicalSolver_Type >                            function_Type;
    typedef BCInterfaceFunctionParser< physicalSolver_Type >                      functionParser_Type;
    typedef OneDFSISolver::solutionPtr_Type                                       solutionPtr_Type;

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
    void updatePhysicalSolverVariables() { std::cout << " !!! WARNING: updatePhysicalSolverVariables() is not defined for the selected solver. !!!" << std::endl; }

    //@}


    //! @name Set Methods
    //@{

    //! Set data for 0D boundary conditions
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void setData( const BCInterfaceData0D& data );

    //! Set data for 1D boundary conditions
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void setData( const BCInterfaceData1D& data );

    //! Set data for 3D boundary conditions
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void setData( const BCInterfaceData3D& data );

    //! Set the physical solver
    /*!
     * @param physicalSolver physical solver
     */
    void setPhysicalSolver( const boost::shared_ptr< PhysicalSolverType >& physicalSolver ) { M_physicalSolver = physicalSolver; }

    //! Set solution
    /*!
     * @param solution The solution container of the 1D problem
     */
    void setSolution( const solutionPtr_Type solution ) { M_solution = solution; }

    //! Set variable function
    /*!
     * @param name name of the variable
     * @param value value of the variable
     */
    void setVariable( const std::string& name, const Real& value ) { functionParser_Type::M_parser->setVariable( name, value ); }

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Create the list of variables for the physical solver.
    /*!
     *  NOTE: A template specialization of this method should be provided for each solver.
     */
    void createAccessList( const BCInterfaceData& /*data*/ ) { std::cout << " !!! WARNING: createAccessList() is not defined for the selected solver. !!!" << std::endl; }

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

    boost::shared_ptr< PhysicalSolverType >    M_physicalSolver;
    solutionPtr_Type                           M_solution;

    OneDFSI::bcSide_Type                       M_side;
    bcFlag_Type                                M_flag;
    std::set< physicalSolverList >             M_list;

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionParserSolver( const BCInterfaceFunctionParserSolver& function );

    BCInterfaceFunctionParserSolver& operator=( const BCInterfaceFunctionParserSolver& function );

    //@}


    //! @name Private Methods
    //@{

    void createFluidMap( std::map< std::string, physicalSolverList >& mapList );
    void createSolidMap( std::map< std::string, physicalSolverList >& mapList );
    void createList( const std::map< std::string, physicalSolverList >& mapList, const BCInterfaceData& data );

    void switchErrorMessage( const std::string& operatorType ) { std::cout << "ERROR: Invalid variable type for " << operatorType << " FunctionSolver" << std::endl; }

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename PhysicalSolverType >
inline BCInterfaceFunctionParser< PhysicalSolverType >* createBCInterfaceFunctionParserSolver()
{
    return new BCInterfaceFunctionParserSolver< PhysicalSolverType > ();
}

// ===================================================
// Constructors
// ===================================================
template< class PhysicalSolverType >
BCInterfaceFunctionParserSolver< PhysicalSolverType >::BCInterfaceFunctionParserSolver() :
        function_Type                    (),
        functionParser_Type              (),
        M_physicalSolver                 (),
        M_solution                       (),
        M_side                           (),
        M_flag                           (),
        M_list                           ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver::BCInterfaceFunctionSolver()" << "\n";
#endif

}

// ===================================================
// Methods
// ===================================================
template< >
inline void
BCInterfaceFunctionParserSolver< OneDFSISolver >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver<FSI>::updatePhysicalSolverVariables" << "\n";
#endif

    // Create/Update variables
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
         // f_ -> FLUID
        case f_timeStep:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_timeStep(): " << M_physicalSolver->physics()->data()->dataTime()->timeStep() << "\n";
#endif
            setVariable( "f_timeStep", M_physicalSolver->physics()->data()->dataTime()->timeStep() );

            break;

        case f_area:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_area(" << static_cast<Real> (M_side) << "): " << M_physicalSolver->boundaryValue( *M_solution, OneDFSI::A, M_side ) << "\n";
#endif
            setVariable( "f_area", M_physicalSolver->boundaryValue( *M_solution, OneDFSI::A, M_side ) );

            break;

        case f_density:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_density: " << M_physicalSolver->physics()->data()->densityRho() << "\n";
#endif
            setVariable( "f_density", M_physicalSolver->physics()->data()->densityRho() );

            break;

        case f_flux:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_flux(" << static_cast<Real> (M_side) << "): " << M_physicalSolver->boundaryValue( *M_solution, OneDFSI::Q, M_side ) << "\n";
#endif

            setVariable( "f_flux", M_physicalSolver->boundaryValue( *M_solution, OneDFSI::Q, M_side ) );

            break;

        case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_pressure(" << static_cast<Real> (M_side) << "): " << M_physicalSolver->boundaryValue( *M_solution, OneDFSI::P, M_side ) << "\n";
#endif

            setVariable( "f_pressure", M_physicalSolver->boundaryValue( *M_solution, OneDFSI::P, M_side ) );

            break;

        case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_viscosity: " << M_physicalSolver->physics()->data()->viscosity() << "\n";
#endif
            setVariable( "f_viscosity", M_physicalSolver->physics()->data()->viscosity() );

            break;

        case f_venousPressure:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_venousPressure: " << M_physicalSolver->physics()->data()->venousPressure() << "\n";
#endif
            setVariable( "f_venousPressure", M_physicalSolver->physics()->data()->venousPressure() );

            break;

        // s_ -> SOLID
        case s_density:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              s_density: " << M_physicalSolver->physics()->data()->densityWall() << "\n";
#endif

            setVariable( "s_density", M_physicalSolver->physics()->data()->densityWall() );

            break;

        case s_poisson:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              s_poisson: " << M_physicalSolver->physics()->data()->poisson() << "\n";
#endif

            setVariable( "s_poisson", M_physicalSolver->physics()->data()->poisson() );

            break;

        case s_thickness:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              s_thickness: " << M_physicalSolver->physics()->data()->thickness( M_physicalSolver->boundaryDOF( M_side ) ) << "\n";
#endif

            setVariable( "s_thickness", M_physicalSolver->physics()->data()->thickness( M_physicalSolver->boundaryDOF( M_side ) ) );

            break;

        case s_young:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              s_young: " << M_physicalSolver->physics()->data()->young() << "\n";
#endif

            setVariable( "s_young", M_physicalSolver->physics()->data()->young() );

            break;

        case s_externalPressure:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              s_externalPressure: " << M_physicalSolver->physics()->data()->externalPressure() << "\n";
#endif

            setVariable( "s_externalPressure", M_physicalSolver->physics()->data()->externalPressure() );

            break;

        default:
            switchErrorMessage( "OneDFSIModel_Solver" );

            break;
        }
}

template< >
inline void
BCInterfaceFunctionParserSolver< FSIOperator >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver<FSIOperator>::updatePhysicalSolverVariables" << "\n";
#endif

    // Create/Update variables
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
        // f_ -> FLUID
        case f_timeStep:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_timeStep(): " << M_physicalSolver->data().dataFluid()->dataTime()->timeStep() << "\n";
#endif
            setVariable( "f_timeStep", M_physicalSolver->data().dataFluid()->dataTime()->timeStep() );

            break;

        case f_area:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_area(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->fluid().area( M_flag ) << "\n";
#endif
            setVariable( "f_area", M_physicalSolver->fluid().area( M_flag ) );

            break;

        case f_density:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_density: " << M_physicalSolver->fluid().density() << "\n";
#endif
            setVariable( "f_density", M_physicalSolver->fluid().density() );

            break;

        case f_flux:

               if ( M_physicalSolver->isFluid() )
            {
#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_flux(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->fluid().flux( M_flag ) << "\n";

#endif
            setVariable( "f_flux", 0.0 );
            }
            else
            {
#ifdef HAVE_LIFEV_DEBUG
                debugStream( 5023 ) << "                                              f_flux(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->fluid().flux( M_flag, *M_physicalSolver->fluid().solution() ) << "\n";
#endif
                setVariable( "f_flux", M_physicalSolver->fluid().flux( M_flag, *M_physicalSolver->fluid().solution() ) );
            }

            break;

        case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_pressure(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->fluid().pressure( M_flag, *M_physicalSolver->fluid().solution() ) << "\n";
#endif

            setVariable( "f_pressure", M_physicalSolver->fluid().pressure( M_flag, *M_physicalSolver->fluid().solution() ) );

            break;

        case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_viscosity: " << M_physicalSolver->fluid().viscosity() << "\n";
#endif
            setVariable( "f_viscosity", M_physicalSolver->fluid().viscosity() );

            break;

        // s_ -> SOLID
        case s_density:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              s_density: " << M_physicalSolver->solid().rho() << "\n";
#endif

            setVariable( "s_density", M_physicalSolver->solid().rho() );

            break;

        case s_poisson:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              s_poisson: " << M_physicalSolver->solid().poisson() << "\n";
#endif

            setVariable( "s_poisson", M_physicalSolver->solid().poisson(1) );

            break;

        case s_thickness:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              s_thickness: " << M_physicalSolver->solid().thickness() << "\n";
#endif

            setVariable( "s_thickness", M_physicalSolver->solid().thickness() );

            break;

        case s_young:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              s_young: " << M_physicalSolver->solid().young() << "\n";
#endif

            setVariable( "s_young", M_physicalSolver->solid().young(1) );

            break;

        case s_externalPressure:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              s_externalPressure: " << M_physicalSolver->solid().data()->externalPressure() << "\n";
#endif

            setVariable( "s_externalPressure", M_physicalSolver->solid().data()->externalPressure() );

            break;

        default:

            switchErrorMessage( "FSIOperator" );

            break;
        }
}

template< >
inline void
BCInterfaceFunctionParserSolver< OseenSolver< RegionMesh< LinearTetra > > >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver<OseenSolver>::updatePhysicalSolverVariables" << "\n";
#endif

    // Create/Update variables
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
        // f_ -> FLUID
        case f_timeStep:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_timeStep(): " << M_physicalSolver->data()->dataTime()->timeStep() << "\n";
#endif
            setVariable( "f_timeStep", M_physicalSolver->data()->dataTime()->timeStep() );

            break;

        case f_area:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_area(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->area( M_flag ) << "\n";
#endif
            setVariable( "f_area", M_physicalSolver->area( M_flag ) );

            break;

        case f_density:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_density: " << M_physicalSolver->density() << "\n";
#endif
            setVariable( "f_density", M_physicalSolver->density() );

            break;

        case f_flux:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_flux(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->flux( M_flag ) << "\n";
#endif

            setVariable( "f_flux", M_physicalSolver->flux( M_flag ) );

            break;

        case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_pressure(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->pressure( M_flag ) << "\n";
#endif

            setVariable( "f_pressure", M_physicalSolver->pressure( M_flag ) );

            break;

        case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_viscosity: " << M_physicalSolver->viscosity() << "\n";
#endif
            setVariable( "f_viscosity", M_physicalSolver->viscosity() );

            break;

        default:

            switchErrorMessage( "OSEEN" );

            break;
        }
}

template< >
inline void
BCInterfaceFunctionParserSolver< OseenSolverShapeDerivative< RegionMesh< LinearTetra > > >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver<OseenSolverShapeDerivative>::updatePhysicalSolverVariables" << "\n";
#endif

    // Create/Update variables
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
        // f_ -> FLUID
        case f_timeStep:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_timeStep(): " << M_physicalSolver->data()->dataTime()->timeStep() << "\n";
#endif
            setVariable( "f_timeStep", M_physicalSolver->data()->dataTime()->timeStep() );

            break;

        case f_area:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_area(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->area( M_flag ) << "\n";
#endif
            setVariable( "f_area", M_physicalSolver->area( M_flag ) );

            break;

        case f_density:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_density(): " << M_physicalSolver->density() << "\n";
#endif
            setVariable( "f_density", M_physicalSolver->density() );

            break;

        case f_flux:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_flux(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->flux( M_flag ) << "\n";
#endif

            setVariable( "f_flux", M_physicalSolver->flux( M_flag ) );

            break;

        case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_pressure(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->pressure( M_flag ) << "\n";
#endif

            setVariable( "f_pressure", M_physicalSolver->pressure( M_flag ) );

            break;

        case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_viscosity(): " << M_physicalSolver->viscosity() << "\n";
#endif
            setVariable( "f_viscosity", M_physicalSolver->viscosity() );

            break;

        default:

            switchErrorMessage( "OSEENSHAPEDERIVATIVE" );

            break;
        }
}

template< >
inline void
#ifdef MULTISCALE_IS_IN_LIFEV
BCInterfaceFunctionParserSolver< Multiscale::MultiscaleData >::updatePhysicalSolverVariables()
#else
BCInterfaceFunctionParserSolver< ZeroDimensionalTemporaryData >::updatePhysicalSolverVariables()
#endif
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver<MultiscaleData>::updatePhysicalSolverVariables" << "\n";
#endif

    // Create/Update variables
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
        // f_ -> FLUID
#ifdef MULTISCALE_IS_IN_LIFEV
        case f_timeStep:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_timeStep(): " << M_physicalSolver->dataTime()->timeStep() << "\n";
#endif
            setVariable( "f_timeStep", M_physicalSolver->dataTime()->timeStep() );

            break;

        case f_density:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_density(): " << M_physicalSolver->fluidDensity() << "\n";
#endif
            setVariable( "f_density", M_physicalSolver->fluidDensity() );

            break;

        case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_viscosity(): " << M_physicalSolver->fluidViscosity() << "\n";
#endif
            setVariable( "f_viscosity", M_physicalSolver->fluidViscosity() );

            break;
#endif
        case f_venousPressure:

#ifdef HAVE_LIFEV_DEBUG
            debugStream( 5023 ) << "                                              f_venousPressure(): " << M_physicalSolver->fluidVenousPressure() << "\n";
#endif
            setVariable( "f_venousPressure", M_physicalSolver->fluidVenousPressure() );

            break;

        default:

            switchErrorMessage( "MultiscaleData" );

            break;
        }
}

// ===================================================
// Set Methods
// ===================================================
template< class PhysicalSolverType >
inline void
BCInterfaceFunctionParserSolver< PhysicalSolverType >::setData( const BCInterfaceData0D& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver::setData( data )" << "\n";
#endif

    M_flag = data.flag();

    functionParser_Type::setData( data );

    createAccessList( data );
}

template< class PhysicalSolverType >
inline void
BCInterfaceFunctionParserSolver< PhysicalSolverType >::setData( const BCInterfaceData1D& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver::setData( data )" << "\n";
#endif

    M_side = data.side();

    functionParser_Type::setData( data );

    createAccessList( data );
}

template< class PhysicalSolverType >
inline void
BCInterfaceFunctionParserSolver< PhysicalSolverType >::setData( const BCInterfaceData3D& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver::setData( data )" << "\n";
#endif

    M_flag = data.flag();

    functionParser_Type::setData( data );

    createAccessList( data );
}

// ===================================================
// Protected Methods
// ===================================================
template< >
inline void
BCInterfaceFunctionParserSolver< OneDFSISolver >::createAccessList( const BCInterfaceData& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver<OneDimensionaSolver>::createAccessList( data )" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createFluidMap( mapList );
    createSolidMap( mapList );
    createList( mapList, data );

    if ( M_physicalSolver.get() )
        updatePhysicalSolverVariables();
}

template< >
inline void
BCInterfaceFunctionParserSolver< FSIOperator >::createAccessList( const BCInterfaceData& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver<FSIOperator>::createAccessList( data )" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createFluidMap( mapList );
    createSolidMap( mapList );
    createList( mapList, data );

    if ( M_physicalSolver.get() )
        updatePhysicalSolverVariables();
}

template< >
inline void
BCInterfaceFunctionParserSolver< OseenSolver< RegionMesh< LinearTetra > > >::createAccessList( const BCInterfaceData& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver<OseenSolver>::createAccessList( data )" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createFluidMap( mapList );
    createList( mapList, data );

    if ( M_physicalSolver.get() )
        updatePhysicalSolverVariables();
}

template< >
inline void
BCInterfaceFunctionParserSolver< OseenSolverShapeDerivative< RegionMesh< LinearTetra > > >::createAccessList( const BCInterfaceData& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver<OseenSolverShapeDerivative>::createAccessList( data )" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createFluidMap( mapList );
    createList( mapList, data );

    if ( M_physicalSolver.get() )
        updatePhysicalSolverVariables();
}

template< >
inline void
#ifdef MULTISCALE_IS_IN_LIFEV
BCInterfaceFunctionParserSolver< Multiscale::MultiscaleData >::createAccessList( const BCInterfaceData& data )
#else
BCInterfaceFunctionParserSolver< ZeroDimensionalTemporaryData >::createAccessList( const BCInterfaceData& data )
#endif
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 5023 ) << "BCInterfaceFunctionSolver<MultiscaleData>::createAccessList( data )" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createFluidMap( mapList );
    createList( mapList, data );

    if ( M_physicalSolver.get() )
        updatePhysicalSolverVariables();
}

// ===================================================
// Private Methods
// ===================================================
template< class PhysicalSolverType >
inline void
BCInterfaceFunctionParserSolver< PhysicalSolverType >::createFluidMap( std::map< std::string, physicalSolverList >& mapList )
{
    mapList["f_timeStep"]       = f_timeStep;
    mapList["f_area"]           = f_area;
    mapList["f_density"]        = f_density;
    mapList["f_flux"]           = f_flux;
    mapList["f_pressure"]       = f_pressure;
    mapList["f_viscosity"]      = f_viscosity;
    mapList["f_venousPressure"] = f_venousPressure;
}

template< class PhysicalSolverType >
inline void
BCInterfaceFunctionParserSolver< PhysicalSolverType >::createSolidMap( std::map< std::string, physicalSolverList >& mapList )
{
    mapList["s_density"]          = s_density;
    mapList["s_poisson"]          = s_poisson;
    mapList["s_thickness"]        = s_thickness;
    mapList["s_young"]            = s_young;
    mapList["s_externalPressure"] = s_externalPressure;
}

template< class PhysicalSolverType >
inline void
BCInterfaceFunctionParserSolver< PhysicalSolverType >::createList( const std::map< std::string, physicalSolverList >& mapList, const BCInterfaceData& data )
{
    M_list.clear();
    for ( typename std::map< std::string, physicalSolverList >::const_iterator j = mapList.begin(); j != mapList.end(); ++j )
        if ( boost::find_first( data.baseString(), j->first ) )
            M_list.insert( j->second );
}

} // Namespace LifeV

#endif /* BCInterfaceFunctionParserSolver_H */
