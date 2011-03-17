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
 *  @brief File containing the BCInterface1D class
 *
 *  @date 10-05-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface1D_H
#define BCInterface1D_H 1

#include <lifemc/lifesolver/BCInterfaceDefinitions.hpp>

#include <lifemc/lifesolver/BCInterfaceData.hpp>
#include <lifemc/lifesolver/BCInterface3DFunction.hpp>
#include <lifemc/lifesolver/BCInterface3DFunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface3DFunctionSolver.hpp>
#include <lifemc/lifesolver/BCInterface3DFunctionFileSolver.hpp>
#include <lifemc/lifesolver/BCInterface1DFunctionDefault.hpp>

namespace LifeV
{

//! BCInterface1D - LifeV Interface to load Boundary Conditions completely from a GetPot file
/*!
 *  @author Cristiano Malossi
 *
 *  This class allows to impose boundary conditions completely from a file.
 *
 *  <b>EXAMPLE - DATA FILE</b>
 *
 *  In the GetPot data file, BCInterface1D reads a new section: [boundary_conditions].
 *
 *  Inside the new section there is a list of condition which correspond to other sub-section
 *  with the same name. The list must be inside the apex ' '.
 *
 *  Each condition has a similar structure; here there is an example:
 *
 *  [InFlow]               <br>
 *  type       = Q         <br>
 *  side       = left      <br>
 *  line       = first     <br>
 *  function   = '3*0.03*(1/4-(x^2+y^2)' <br>
 *
 *  NOTE: All the parameters are case sensitive.
 *
 *  type - can be: A, Q, W1, W2, P
 *  side - can be: left, right
 *  line - can be: first, second.
 *  function - contains the function. See BCInterfaceFunction1D for more details about the syntax.
 *
 *  <b>NOTE:</b>
 *
 *  The string "function" represent the base module and can be replaced by other expanded modules.
 *  Up to now we have the following modules for function:
 *
 *  - function
 *  - functionFile
 *  - OperatorFunction
 *  - OperatorFunctionFile
 *  - Default_1D
 *
 *  To see some example look at test_fsi.
 *
 *  <b>EXAMPLE - HOW TO USE</b>
 *
 *  Here there is a short example on how to use it.
 *
 *  1) You can define your BCInterface1D class in a shared pointer:
 *     boost::shared_ptr<BCInterface1D> 	M_fluidBC;
 *
 *  2) Build the BCInterface1D using empty constructor;
 *
 *  3) If you have operator conditions you have to give the operator to access variables
 *     M_fluidBC->setPhysicalSolver( M_fsi->FSIOper() );
 *
 *  4) Then you can fill the handler from a file and a section (this can be done for multiple files & sections)
 *     M_fluidBC->fillHandler( "fileName.dat", "fluid" );
 *
 *  5) Finally, to get the handler you can use:
 *     M_fluidBC->handler();
 *
 *  NOTE:
 *
 *  a) You can add manually more conditions by using setBC() after the call to buildHandler() function.
 *     In this case you have to manually set the TOTAL number of boundary conditions
 *     by using setHandlerParameters() function BEFORE building the handler.
 */
template< class PhysicalSolverType = OneDimensionalSolver >
class BCInterface1D
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                                        physicalSolver_Type;
    typedef baseList1D_Type                                                                           bcBaseList_Type;
    typedef BCInterfaceData                                                                           data_Type;

    typedef FactorySingleton< Factory< BCInterface3DFunction< physicalSolver_Type > , bcBaseList_Type > >    factoryFunction_Type;

    typedef OneDimensionalBCHandler                                                                   bcHandler_Type;
    typedef boost::shared_ptr< bcHandler_Type >                                                       bcHandlerPtr_Type;

    typedef bcHandler_Type::solutionPtr_Type                                                          solutionPtr_Type;
    typedef bcHandler_Type::fluxPtr_Type                                                              fluxPtr_Type;
    typedef bcHandler_Type::sourcePtr_Type                                                            sourcePtr_Type;

    typedef bcHandler_Type::vectorPtrContainer_Type                                                   vectorPtrContainer_Type;

    typedef std::vector< boost::shared_ptr< BCInterface3DFunction< physicalSolver_Type > > >          vectorFunction_Type;
    typedef std::vector< boost::shared_ptr< BCInterface1DFunctionDefault< physicalSolver_Type > > >   vectorDefaultFunction_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface1D();

    //! Destructor
    virtual ~BCInterface1D() {}

    //@}


    //! @name Methods
    //@{

    //! Create the bcHandler.
    void createHandler() { M_handler.reset( new bcHandler_Type() ); }

    //! Fill the bcHandler with the BC provided in the file.
    /*!
     * @param fileName Name of the data file
     * @param dataSection Subsection inside [boundary_conditions]
     */
    void fillHandler( const std::string& fileName, const std::string& dataSection );

    //! Read a specific boundary condition from a file and add it to the data container
    /*!
     * @param fileName Name of the data file
     * @param dataSection section in the data file
     * @param name name of the boundary condition
     */
    void readBC( const std::string& fileName, const std::string& dataSection, const bcName_Type& name )
    {
        M_data.readBC1D( fileName, dataSection, name );
    }

    //! Insert the current boundary condition in the BChandler
    void insertBC();

    //! Update the variables inside the physical solver
    void updatePhysicalSolverVariables();

    //@}


    //! @name External interface for BCHandler functions
    //@{

    //! Add a Boundary Condition using the standard interface of the BCHandler
    /*!
     * @param side side of the condition
     * @param line line of the condition
     * @param type type of the condition
     * @param base base of the condition
     */
    template< class BCBaseType >
    void setBC( const OneDimensional::bcSide_Type& bcSide, const OneDimensional::bcLine_Type& bcLine, const OneDimensional::bcType_Type& bcType, const BCBaseType& base )
    {
        M_handler->setBC( bcSide, bcLine, bcType, base );
    }

    //@}


    //! @name Set Methods
    //@{

    //! Set a physical solver
    /*!
     * @param physicalSolver physical solver
     */
    void setPhysicalSolver( const boost::shared_ptr< physicalSolver_Type >& physicalSolver );

    //! Set an Handler
    /*!
     * @param handler BCHandler
     */
    void setHandler( const bcHandlerPtr_Type& handler ) { M_handler = handler; }

    //! Set the solution for the members that need it
    /*!
     * @param flux flux
     * @param source source
     */
    void setFluxSource( const fluxPtr_Type& flux, const sourcePtr_Type& source );

    //! Set the solution for the members that need it
    /*!
     * @param solution solution
     */
    void setSolution( const solutionPtr_Type& solution );

    //@}


    //! @name Get Methods
    //@{

    //! Get the shared_ptr to the BCHandler
    /*!
     * @return the pointer to the BCHandler
     */
    const bcHandlerPtr_Type& handler() { return M_handler; }

    //! Get the data container
    /*!
     * @return the data container
     */
    data_Type& dataContainer() { return M_data; }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface1D( const BCInterface1D& interface1D );

    BCInterface1D& operator=( const BCInterface1D& interface1D );

    //@}


    //! @name Private Methods
    //@{

    template< class BCInterfaceBaseType >
    void createFunction( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector );

    template< class BCInterfaceBaseType >
    void createFunction( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector, const bcBaseList_Type& physicalSolver );

    template< class BCBaseType >
    void addBcToHandler( BCBaseType& base );

    //@}

    // Handler and parameters
    bcHandlerPtr_Type                        M_handler;

    // Data
    data_Type                                M_data;

    // Functions
    vectorFunction_Type                      M_vectorFunction;

    // Default Functions
    vectorDefaultFunction_Type               M_vectorDefaultFunction1D;
};

// ===================================================
// Constructors & Destructor
// ===================================================
template< class PhysicalSolverType >
BCInterface1D< PhysicalSolverType >::BCInterface1D( ) :
        M_handler                 (),
        M_data                    (),
        M_vectorFunction          (),
        M_vectorDefaultFunction1D ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface1D::BCInterface1D" << "\n";
#endif

    //Factory registration
    factoryFunction_Type::instance().registerProduct( BCI1DFunction,           &createBCInterface3DFunction< physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct( BCI1DFunctionFile,       &createBCInterface3DFunctionFile< physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct( BCI1DFunctionSolver,     &createBCInterface3DFunctionSolver< physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct( BCI1DFunctionFileSolver, &createBCInterface3DFunctionFileSolver< physicalSolver_Type > );
}

// ===================================================
// Methods
// ===================================================
template< class PhysicalSolverType >
void
BCInterface1D< PhysicalSolverType >::fillHandler( const std::string& fileName, const std::string& dataSection )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface1D::fillHandler\n";
#endif

    GetPot dataFile( fileName );
    for ( UInt i( 0 ); i < dataFile.vector_variable_size( ( dataSection + "/boundary_conditions/list" ).c_str() ); ++i )
    {
        M_data.readBC1D( fileName,
                         dataSection + "/boundary_conditions/",
                         dataFile( ( dataSection + "/boundary_conditions/list" ).c_str(), " ", i )
                       );

        insertBC();
    }
}

template< class PhysicalSolverType >
inline void
BCInterface1D< PhysicalSolverType >::insertBC()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface1D::insertBC\n";
#endif

    switch ( M_data.base1D().second )
    {
    case BCI1DFunction:
    case BCI1DFunctionFile:
    case BCI1DFunctionSolver:
    case BCI1DFunctionFileSolver:
    {
        createFunction( M_vectorFunction, M_data.base1D().second );

        OneDimensionalBCFunction base;
        M_vectorFunction.back()->assignFunction( base );

        addBcToHandler( base );

        break;
    }
    case BCI1DFunctionDefault:
    {
        createFunction( M_vectorDefaultFunction1D );

        OneDimensionalBCFunction base;
        M_vectorDefaultFunction1D.back()->assignFunction( base );

        addBcToHandler( base );

        break;
    }
    }
}

template< class PhysicalSolverType >
void
BCInterface1D< PhysicalSolverType >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface1D::updatePhysicalSolverVariables\n";
#endif

    for ( UInt i( 0 ); i < M_vectorFunction.size(); ++i )
    {
        BCInterface3DFunctionSolver< physicalSolver_Type > *physicalSolver =
            dynamic_cast< BCInterface3DFunctionSolver< physicalSolver_Type > * > ( &( *M_vectorFunction[i] ) );

        if ( physicalSolver != 0 )
            physicalSolver->updatePhysicalSolverVariables();
    }
}

// ===================================================
// Set Methods
// ===================================================
template< class PhysicalSolverType >
void
BCInterface1D< PhysicalSolverType >::setPhysicalSolver( const boost::shared_ptr< physicalSolver_Type >& physicalSolver )
{
    //for ( typename vectorFunction_Type::const_iterator i = M_vectorFunction.begin() ; i < M_vectorFunction.end() ; ++i )
    for ( UInt i( 0 ); i < M_vectorFunction.size(); ++i )
    {
        BCInterface3DFunctionSolver< physicalSolver_Type > *castedOperator =
            dynamic_cast < BCInterface3DFunctionSolver< physicalSolver_Type >* > ( &( *M_vectorFunction[i] ) );

        if ( castedOperator != 0 )
            castedOperator->setPhysicalSolver( physicalSolver );
    }
}

template< class PhysicalSolverType >
void
BCInterface1D< PhysicalSolverType >::setFluxSource( const fluxPtr_Type& flux, const sourcePtr_Type& source )
{
    for ( typename vectorDefaultFunction_Type::const_iterator i = M_vectorDefaultFunction1D.begin() ; i < M_vectorDefaultFunction1D.end() ; ++i )
        ( *i )->setFluxSource( flux, source );

    M_handler->setFluxSource( flux, source );
}

template< class PhysicalSolverType >
void
BCInterface1D< PhysicalSolverType >::setSolution( const solutionPtr_Type& solution )
{
    //for ( typename vectorFunction_Type::const_iterator i = M_vectorFunction.begin() ; i < M_vectorFunction.end() ; ++i )
    for ( UInt i( 0 ); i < M_vectorFunction.size(); ++i )
    {
        BCInterface3DFunctionSolver< physicalSolver_Type > *castedOperator =
            dynamic_cast < BCInterface3DFunctionSolver< physicalSolver_Type >* > ( &( *M_vectorFunction[i] ) );

        if ( castedOperator != 0 )
            castedOperator->setSolution( solution );
    }

    for ( typename vectorDefaultFunction_Type::const_iterator i = M_vectorDefaultFunction1D.begin() ; i < M_vectorDefaultFunction1D.end() ; ++i )
        ( *i )->setSolution( solution );

    M_handler->setSolution( solution );
}

// ===================================================
// Private Methods
// ===================================================
template< class PhysicalSolverType > template< class BCInterfaceBaseType >
inline void
BCInterface1D< PhysicalSolverType >::createFunction( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector )
{
    boost::shared_ptr< BCInterfaceBaseType > function( new BCInterfaceBaseType( M_data ) );
    baseVector.push_back( function );
}

template< class PhysicalSolverType > template< class BCInterfaceBaseType >
inline void
BCInterface1D< PhysicalSolverType >::createFunction( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector, const bcBaseList_Type& physicalSolver )
{
    boost::shared_ptr< BCInterfaceBaseType > function( factoryFunction_Type::instance().createObject( physicalSolver, M_data.mapBase1D() ) );

    function->setData( M_data );

    baseVector.push_back( function );
}

template< class PhysicalSolverType > template< class BCBaseType >
inline void
BCInterface1D< PhysicalSolverType >::addBcToHandler( BCBaseType& base )
{
    if ( !M_handler.get() ) // If BCHandler has not been created yet, we do it now
        createHandler();

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface1D::addBCManager" << "\n\n";
#endif

    M_handler->setBC( M_data.side(), M_data.line(), M_data.quantity(), base );
}

} // Namespace LifeV

#endif /* BCInterface1D_H */
