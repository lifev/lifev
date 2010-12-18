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
 *  @brief File containing the BCInterface3D class
 *
 *  @date 01-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface3D_H
#define BCInterface3D_H 1

#include <lifemc/lifesolver/BCInterface3DDefinitions.hpp>

#include <lifemc/lifesolver/BCInterface3DData.hpp>
#include <lifemc/lifesolver/BCInterface3DFunction.hpp>
#include <lifemc/lifesolver/BCInterface3DFunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface3DFunctionSolver.hpp>
#include <lifemc/lifesolver/BCInterface3DFunctionFileSolver.hpp>
#include <lifemc/lifesolver/BCInterface3DFunctionFSI.hpp>

namespace LifeV
{

//! BCInterface3D - LifeV Interface to load Boundary Conditions completely from a GetPot file
/*!
 *  @author Cristiano Malossi
 *
 *  This class allows to impose boundary conditions completely from a file.
 *
 *  <b>EXAMPLE - DATA FILE</b>
 *
 *  In the GetPot data file, inside each subsection ([fluid], [solid], etc...) BCInterface3D reads a new section: [boundary_conditions].
 *
 *  Inside the new section there is a list of condition which correspond to other sub-section
 *  with the same name. The list must be inside the apex ' '.
 *
 *  Each condition has a similar structure; here there is an example:
 *
 *  [InFlow]               <br>
 *  type       = Essential <br>
 *  flag       = 2         <br>
 *  mode       = Full      <br>
 *  component  = 3         <br>
 *  function   = '[0, 0, 3*0.03*(1/4-(x^2+y^2)]' <br>
 *
 *  NOTE: All the parameters are case sensitive.
 *
 *  type - can be: Essential Natural Mixte Flux
 *  flag - contains the flag
 *  mode - can be: Full Component Scalar Tangential Normal.
 *  component - if mode is Scalar, Tangential or Normal it is missing.
 *              if mode is Component it contains the ID of the component (or of the components list inside apex)
 *              if mode is Full it contains the total number of components
 *  function - contains the function. See BCInterfaceFunction for more details about the syntax.
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
 *  - FSI
 *
 *  To see some example look at test_fsi.
 *
 *  <b>EXAMPLE - HOW TO USE</b>
 *
 *  Here there is a short example on how to use it.
 *
 *  1) You can define your BCInterface class in a shared pointer:
 *     boost::shared_ptr<BCInterface> 	M_fluidBC;
 *
 *  2) Build the BCInterface using empty constructor;
 *
 *  3) If you have operator conditions you have to give the operator to access variables
 *     M_fluidBC->setPhysicalSolver( M_fsi->FSIOper() );
 *
 *  4) Then you can fill the handler from a file and a section (this can be done for multiple files & sections)
 *     M_fluidBC->fillHandler( "fileName.dat", "fluid" );
 *
 *  5) Finally, to get the handler you can use:
 *     M_fluidBC->handler();
 */

template< class PhysicalSolverType >
class BCInterface3D
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                                        physicalSolver_Type;
    typedef baseList3D_Type                                                                           bcBaseList_Type;
    typedef BCInterface3DData                                                                         data_Type;

    typedef singleton< factory< BCInterface3DFunction< physicalSolver_Type > , bcBaseList_Type > >    factoryFunction_Type;

    typedef BCHandler                                                                                 bcHandler_Type;
    typedef boost::shared_ptr< bcHandler_Type >                                                       bcHandlerPtr_Type;

    typedef std::vector< boost::shared_ptr< BCInterface3DFunction< physicalSolver_Type > > >          vectorFunction_Type;
    typedef std::vector< boost::shared_ptr< BCFunctionDirectional > >                                 vectorFunctionDirectional_Type;
    typedef std::vector< boost::shared_ptr< BCInterface3DFunctionFSI< physicalSolver_Type > > >       vectorFSI_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface3D();

    //! Destructor
    virtual ~BCInterface3D() {}

    //@}


    //! @name Methods
    //@{

    //! Update the variables inside the physical solver
    void updatePhysicalSolverVariables();

    //! Create the bcHandler.
    void createHandler() { M_handler.reset( new bcHandler_Type( ) ); }

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
    void readBC( const std::string& fileName, const std::string& dataSection, const bcName_Type& name ) { M_data.readBC( fileName, dataSection, name ); }

    //! Insert the current boundary condition in the BChandler
    void insertBC() { buildBase(); }

    //@}


    //! @name External interface for BCHandler functions
    //@{

    //! Add a Boundary Condition using the standard interface of the BCHandler
    /*!
     * @param name name of the condition
     * @param flag list of flags
     * @param type type of the condition
     * @param mode mode of the condition
     * @param base base of the condition
     */
    template< class BCBaseType >
    void addBC( const bcName_Type& name, const bcFlag_Type& flag, const bcType_Type& type, const bcMode_Type& mode, BCBaseType& base );

    //! Add a Boundary Condition with component using the standard interface of the BCHandler
    /*!
     * @param name name of the condition
     * @param flag list of flags
     * @param type type of the condition
     * @param mode mode of the condition
     * @param base base of the condition
     * @param comp component of the condition
     */
    template< class BCBaseType, class BCCompType >
    void addBC( const bcName_Type& name, const bcFlag_Type& flag, const bcType_Type& type, const bcMode_Type& mode, BCBaseType& base, const BCCompType& comp );

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

    BCInterface3D( const BCInterface3D& bcInterface3D );

    BCInterface3D& operator=( const BCInterface3D& bcInterface3D );

    //@}


    //! @name Private Methods
    //@{

    void buildBase();

    template< class BCInterfaceBaseType >
    void addBase( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector );

    template< class BCInterfaceBaseType >
    void addBase( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector, const bcBaseList_Type& physicalSolver );

    // This method should be removed: it is a workaround due to legacy of LifeV BC.
    void addBCManager( BCVectorInterface& base );

    template< class BCBaseType >
    void addBCManager( BCBaseType& base );

    //@}

    // Handler and parameters
    bcHandlerPtr_Type               M_handler;

    // Data
    data_Type                       M_data;

    // Functions
    vectorFunction_Type             M_vectorFunction;

    // Functions Directions
    vectorFunctionDirectional_Type  M_vectorFunctionDirection;

    // FSI Functions
    vectorFSI_Type                  M_vectorFSI;
};

// ===================================================
// Constructors & Destructor
// ===================================================
template< class PhysicalSolverType >
BCInterface3D< PhysicalSolverType >::BCInterface3D( ) :
        M_handler                 (),
        M_data                    (),
        M_vectorFunction          (),
        M_vectorFunctionDirection (),
        M_vectorFSI               ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface3D::BCInterface3D------------------------------" << "\n";
#endif

    //Factory registration
    factoryFunction_Type::instance().registerProduct( BCI3DFunction,           &createBCInterface3DFunction< physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct( BCI3DFunctionFile,       &createBCInterface3DFunctionFile< physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct( BCI3DFunctionSolver,     &createBCInterface3DFunctionSolver< physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct( BCI3DFunctionFileSolver, &createBCInterface3DFunctionFileSolver< physicalSolver_Type > );
}

// ===================================================
// Methods
// ===================================================
template< class PhysicalSolverType >
void
BCInterface3D< PhysicalSolverType >::fillHandler( const std::string& fileName, const std::string& dataSection )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface3D::buildHandler\n";
#endif

    GetPot DataFile( fileName );
    for ( UInt i( 0 ); i < DataFile.vector_variable_size( ( dataSection + "/boundary_conditions/list" ).c_str() ); ++i )
    {
        M_data.readBC( fileName,
                       dataSection + "/boundary_conditions/",
                       DataFile( ( dataSection + "/boundary_conditions/list" ).c_str(), " ", i )
                     );

        buildBase();
    }
}

template< class PhysicalSolverType >
void
BCInterface3D< PhysicalSolverType >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface3D::UpdateOperatorVariables\n";
#endif

    for ( UInt i( 0 ); i < M_vectorFunction.size(); ++i )
    {
        BCInterface3DFunctionSolver< physicalSolver_Type > *physicalSolver =
            dynamic_cast< BCInterface3DFunctionSolver< physicalSolver_Type > * > ( &( *M_vectorFunction[i] ) );

        if ( physicalSolver != 0 )
            physicalSolver->updatePhysicalSolverVariables();
    }
}

template< class PhysicalSolverType > template< class BCBaseType >
void
BCInterface3D< PhysicalSolverType >::addBC( const bcName_Type& name, const bcFlag_Type& flag, const bcType_Type& type, const bcMode_Type& mode, BCBaseType& base )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface3D::addBC (without component)" << "\n\n";
#endif

    M_handler->addBC( name, flag, type, mode, base );
}

template< class PhysicalSolverType > template< class BCBaseType, class BCCompType >
void
BCInterface3D< PhysicalSolverType >::addBC( const bcName_Type& name, const bcFlag_Type& flag, const bcType_Type& type, const bcMode_Type& mode, BCBaseType& base, const BCCompType& comp )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface3D::addBC (with component)" << "\n\n";
#endif

    M_handler->addBC( name, flag, type, mode, base, comp );
}

// ===================================================
// Set Methods
// ===================================================
template< class PhysicalSolverType >
void
BCInterface3D< PhysicalSolverType >::setPhysicalSolver( const boost::shared_ptr< physicalSolver_Type >& physicalSolver )
{
    //for ( typename vectorFunction_Type::const_iterator i = M_vectorFunction.begin() ; i < M_vectorFunction.end() ; ++i )
    for ( UInt i( 0 ); i < M_vectorFunction.size(); ++i )
    {
        BCInterface3DFunctionSolver< physicalSolver_Type > *castedOperator =
            dynamic_cast < BCInterface3DFunctionSolver< physicalSolver_Type >* > ( &( *M_vectorFunction[i] ) );

        if ( castedOperator != 0 )
            castedOperator->setPhysicalSolver( physicalSolver );
    }

    for ( typename vectorFSI_Type::const_iterator i = M_vectorFSI.begin() ; i < M_vectorFSI.end() ; ++i )
    {
        ( *i )->checkMethod( physicalSolver );
        ( *i )->exportData( M_data );
        addBCManager( ( *i )->base() );
    }
}


// ===================================================
// Private Methods
// ===================================================
template< class PhysicalSolverType >
inline void
BCInterface3D< PhysicalSolverType >::buildBase()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface3D::BuildBase\n";
#endif

    switch ( M_data.base().second )
    {
    case BCI3DFunction:
    case BCI3DFunctionFile:
    case BCI3DFunctionSolver:
    case BCI3DFunctionFileSolver:

        addBase( M_vectorFunction, M_data.base().second );

        addBCManager( M_vectorFunction.back()->base() );

        break;

    case BCI3DFunctionFSI:

        addBase( M_vectorFSI );

        break;
    }
}

template< class PhysicalSolverType > template< class BCInterfaceBaseType >
inline void
BCInterface3D< PhysicalSolverType >::addBase( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector )
{
    boost::shared_ptr< BCInterfaceBaseType > Function( new BCInterfaceBaseType( M_data ) );
    baseVector.push_back( Function );
}

template< class PhysicalSolverType > template< class BCInterfaceBaseType >
inline void
BCInterface3D< PhysicalSolverType >::addBase( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector, const bcBaseList_Type& physicalSolver )
{
    boost::shared_ptr< BCInterfaceBaseType > function( factoryFunction_Type::instance().createObject( physicalSolver ) );

    //!\todo pass a std::string to the factories
    //boost::shared_ptr< BCInterfaceBaseType > function( factoryFunction_Type::instance().createObject( "physicalSolver" ) );

    function->setData( M_data );

    baseVector.push_back( function );
}

template< class PhysicalSolverType >
inline void
BCInterface3D< PhysicalSolverType >::addBCManager( BCVectorInterface& base )
{
    if ( !M_handler.get() ) // If BCHandler has not been created yet, we do it now
        createHandler();

    switch ( M_data.mode() )
    {
    case Scalar:
    case Normal:
    case Tangential:
    case Directional:
    case Component:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::AddBCManager                              others" << "\n\n";
#endif

        std::cout << "ERROR: Scalar, Normal, Tangential, Directional, Component NOT AVAILABLE FOR FSI BC" << std::endl;
        break;

    case Full:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::AddBCManager                              Full" << "\n\n";
#endif

        M_handler->addBC( M_data.name(), M_data.flag(), M_data.type(), M_data.mode(), base, M_data.comN() );

        break;
    }
}

template< class PhysicalSolverType > template< class BCBaseType >
inline void
BCInterface3D< PhysicalSolverType >::addBCManager( BCBaseType& base )
{
    if ( !M_handler.get() ) // If BCHandler has not been created yet, we do it now
        createHandler();

    switch ( M_data.mode() )
    {
    case Scalar:
    case Normal:
    case Tangential:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::AddBCManager                              Scalar, Normal, Tangential" << "\n\n";
#endif

        M_handler->addBC( M_data.name(), M_data.flag(), M_data.type(), M_data.mode(), base );

        break;

    case Directional:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::AddBCManager                              Directional" << "\n\n";
#endif
        {
            // Parameters for direction BC
            M_data.setName( M_data.name() + "_direction" );
            M_data.setBase( make_pair( "function", BCI3DFunction ) );
            M_data.setBaseString( M_data.direction() );

            // Directional field
            addBase( M_vectorFunction, M_data.base().second );

            // Directional base
            boost::shared_ptr< BCFunctionDirectional > directionalBase( new BCFunctionDirectional( base.Function(), M_vectorFunction.back()->base().Function() ) );
            M_vectorFunctionDirection.push_back( directionalBase );

            M_handler->addBC( M_data.name(), M_data.flag(), M_data.type(), M_data.mode(), *M_vectorFunctionDirection.back() );
        }

        break;

    case Full:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::AddBCManager                              Full" << "\n\n";
#endif

        M_handler->addBC( M_data.name(), M_data.flag(), M_data.type(), M_data.mode(), base, M_data.comN() );

        break;

    case Component:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::AddBCManager                              Component" << "\n\n";
#endif

        M_handler->addBC( M_data.name(), M_data.flag(), M_data.type(), M_data.mode(), base, M_data.comV() );

        break;
    }
}

} // Namespace LifeV

#endif /* BCInterface3D_H */
