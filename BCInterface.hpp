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
 *  @brief BCInterface
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 01-04-2009
 */

#ifndef BCInterface_H
#define BCInterface_H 1

#include <lifemc/lifesolver/BCInterface_Definitions.hpp>

#include <lifemc/lifesolver/BCInterface_Data.hpp>
#include <lifemc/lifesolver/BCInterface_Function.hpp>
#include <lifemc/lifesolver/BCInterface_FunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface_OperatorFunction.hpp>
#include <lifemc/lifesolver/BCInterface_OperatorFunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface_FSI.hpp>

namespace LifeV {

//! BCInterface - LifeV Interface to load Boundary Conditions completely from a GetPot file
/*!
 *  @author Cristiano Malossi
 *
 *  This class allows to impose boundary conditions completely from a file.
 *
 *  <b>EXAMPLE - DATA FILE</b>
 *
 *  In the GetPot data file, inside each subsection ([fluid], [solid], etc...) BCInterface reads a new section: [boundary_conditions].
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
 *     M_fluidBC->setOperator( M_fsi->FSIOper() );
 *
 *  4) Then you can fill the handler from a file and a section (this can be done for multiple files & sections)
 *     M_fluidBC->FillHandler( "FileName.dat", "fluid" );
 *
 *  5) Finally, to get the handler you can use:
 *     M_fluidBC->GetHandler();
 *
 *  NOTE:
 *
 *  a) You can add manually more conditions by using addBC() after the call to buildHandler() function.
 *     In this case you have to manually set the TOTAL number of boundary conditions
 *     by using setHandlerParameters() function BEFORE building the handler.
 */
template< class Operator >
class BCInterface
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface_BaseList                                                       BCBaseList_Type;

    typedef singleton< factory< BCInterface_Function< Operator > , BCBaseList_Type > > FactoryBCInterface_Function;

    typedef BCHandler                                                                  BCHandler_Type;
    typedef boost::shared_ptr< BCHandler_Type >                                        BCHandler_PtrType;

    typedef BCInterface_Data< Operator >                                               Data_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface();

    //! Copy constructor
    /*!
     * @param interface	- BCInterface
     */
    BCInterface( const BCInterface& interface );

    //! Destructor
    ~BCInterface() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param interface BCInterface
     * @return reference to a copy of the class
     */
    BCInterface& operator=( const BCInterface& interface );

    //@}


    //! @name Methods
    //@{

    //! Update the variables inside the operator
    void UpdateOperatorVariables();

    //! Create the bcHandler.
    void CreateHandler();

    //! Fill the bcHandler with the BC provided in the file.
    /*!
     * @param FileName Name of the data file
     * @param dataSection Subsection inside [boundary_conditions]
     */
    void FillHandler( const std::string& FileName, const std::string& dataSection );

    //! Read a specific boundary condition from a file and add it to the data container
    /*!
     * @param FileName Name of the data file
     * @param dataSection section in the data file
     * @param name name of the boundary condition
     */
    void ReadBC( const std::string& FileName,
                 const std::string& dataSection,
                 const BCName&      name );

    //! Insert the current boundary condition in the BChandler
    void InsertBC();

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
    template< class BCBase >
    void addBC( const BCName& name,
                const BCFlag& flag,
                const BCType& type,
                const BCMode& mode,
                      BCBase& base );

    //! Add a Boundary Condition with component using the standard interface of the BCHandler
    /*!
     * @param name name of the condition
     * @param flag list of flags
     * @param type type of the condition
     * @param mode mode of the condition
     * @param base base of the condition
     * @param comp component of the condition
     */
    template< class BCBase, class BCComp >
    void addBC( const BCName& name,
                const BCFlag& flag,
                const BCType& type,
                const BCMode& mode,
                      BCBase& base,
                const BCComp& comp );

    //@}


    //! @name Set Methods
    //@{

    //! Set an operator
    /*!
     * @param Oper operator
     */
    void SetOperator( const boost::shared_ptr< Operator >& Oper );

    //! Set an Handler
    /*!
     * @param handler BCHandler
     */
    void SetHandler( const BCHandler_PtrType& handler );

    //! Set manually Handler parameters: you need it only if you are adding manually some parameters by calling addBC
    /*!
     * @param bcNumber total number of the boundary conditions (files + added manually)
     * @param hint hint
     */
    void SetHandlerParameters( const ID& bcNumber,
                               const BCHandler_Type::BCHints& hint = BCHandler_Type::HINT_BC_NONE );

    //@}


    //! @name Get Methods
    //@{

    //! Get the shared_ptr to the BCHandler
    /*!
     * @return the pointer to the BCHandler
     */
    const BCHandler_PtrType& GetHandler();

    //! Get the data container
    /*!
     * @return the data container
     */
    Data_Type& GetDataContainer();

    //@}

private:

    //! @name Private Methods
    //@{

    inline void BuildBase();

    template< class BCInterfaceBase >
    inline void AddBase( std::vector< boost::shared_ptr< BCInterfaceBase > >& baseVector );

    template< class BCInterfaceBase >
    inline void AddBase(       std::vector< boost::shared_ptr< BCInterfaceBase > >& baseVector,
                         const BCBaseList_Type& Oper );

    // This method should be removed: it is a workaround due to legacy of LifeV BC.
    inline void AddBCManager( BCVectorInterface& base );

    template< class BCBase >
    inline void AddBCManager( BCBase& base );

    //@}

    // Handler and parameters
    ID                              M_bcNumber;
    BCHandler_Type::BCHints         M_hint;
    BCHandler_PtrType               M_handler;

    // Data
    Data_Type                       M_data;

    // Functions
    std::vector< boost::shared_ptr< BCInterface_Function< Operator > > > M_vectorFunction;

    // Functions Directions
    std::vector< boost::shared_ptr< BCFunctionDirectional > >            M_vectorFunctionDirection;

    // FSI Functions
    std::vector< boost::shared_ptr< BCInterface_FSI< Operator > > >      M_vectorFSI;
};

// ===================================================
// Constructors & Destructor
// ===================================================
template< class Operator >
BCInterface< Operator >::BCInterface( ) :
    M_bcNumber                ( 0 ),
    M_hint                    ( BCHandler_Type::HINT_BC_NONE ),
    M_handler                 (),
    M_data                    (),
    M_vectorFunction          (),
    M_vectorFunctionDirection (),
    M_vectorFSI               ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::BCInterface------------------------------" << "\n";
#endif

    //Factory registration
    FactoryBCInterface_Function::instance().registerProduct( BCInterface_function,         &BCInterface_CreateFunction< Operator > );
    FactoryBCInterface_Function::instance().registerProduct( BCInterface_functionFile,     &BCInterface_CreateFunctionFile< Operator > );
    FactoryBCInterface_Function::instance().registerProduct( BCInterface_OPERfunction,     &BCInterface_CreateOperatorFunction< Operator > );
    FactoryBCInterface_Function::instance().registerProduct( BCInterface_OPERfunctionFile, &BCInterface_CreateOperatorFunctionFile< Operator > );
}

template< class Operator >
BCInterface< Operator >::BCInterface( const BCInterface& interface ) :
    M_bcNumber                ( interface.M_bcNumber ),
    M_hint                    ( interface.M_hint ),
    M_handler                 ( interface.M_handler ),
    M_data                    ( interface.M_data ),
    M_vectorFunction          ( interface.M_vectorFunction ),
    M_vectorFunctionDirection ( interface.M_vectorFunctionDirection ),
    M_vectorFSI               ( interface.M_vectorFSI )
{
}

// ===================================================
// Operators
// ===================================================
template< class Operator >
BCInterface< Operator >&
BCInterface< Operator >::operator=( const BCInterface& interface )
{
    if ( this != &interface )
    {
        M_bcNumber                = interface.M_bcNumber;
        M_hint                    = interface.M_hint;
        M_handler                 = interface.M_handler;
        M_data                    = interface.M_data;
        M_vectorFunction          = interface.M_vectorFunction;
        M_vectorFunctionDirection = interface.M_vectorFunctionDirection;
        M_vectorFSI               = interface.M_vectorFSI;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
template< class Operator >
void
BCInterface< Operator >::CreateHandler()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::CreateHandler\n";
#endif

    M_handler.reset( new BCHandler_Type( M_bcNumber, M_hint ) );
}

template< class Operator >
void
BCInterface< Operator >::FillHandler( const std::string& FileName,
                                      const std::string& dataSection )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::buildHandler\n";
#endif

    GetPot DataFile( FileName );
    for ( UInt i( 0 ); i < DataFile.vector_variable_size( ( dataSection + "/boundary_conditions/list" ).c_str() ); ++i )
    {
        M_data.ReadBC( FileName,
                       dataSection + "/boundary_conditions/",
                       DataFile( ( dataSection + "/boundary_conditions/list" ).c_str(), " ", i )
                     );

        BuildBase();
    }
}

template< class Operator >
void
BCInterface< Operator >::ReadBC( const std::string& FileName,
                                 const std::string& dataSection,
                                 const BCName&      name )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::ReadBC\n";
#endif

    M_data.ReadBC( FileName, dataSection, name );
}

template< class Operator >
void
BCInterface< Operator >::InsertBC()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::InsertBC\n";
#endif

    BuildBase();
}

template< class Operator >
void
BCInterface< Operator >::UpdateOperatorVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::UpdateOperatorVariables\n";
#endif

    for ( UInt i( 0 ); i < M_vectorFunction.size(); ++i )
    {
        BCInterface_OperatorFunction< Operator > *Oper =
                dynamic_cast< BCInterface_OperatorFunction< Operator > * > ( &( *M_vectorFunction[i] ) );

        if ( Oper != 0 )
            Oper->UpdateOperatorVariables();
    }
}

template< class Operator > template< class BCBase >
void
BCInterface< Operator >::addBC( const BCName& name,
                                const BCFlag& flag,
                                const BCType& type,
                                const BCMode& mode,
                                      BCBase& base )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::addBC (without component)" << "\n\n";
#endif

    M_handler->addBC( name, flag, type, mode, base );
}

template< class Operator > template< class BCBase, class BCComp >
void
BCInterface< Operator >::addBC( const BCName& name,
                                const BCFlag& flag,
                                const BCType& type,
                                const BCMode& mode,
                                      BCBase& base,
                                const BCComp& comp )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::addBC (with component)" << "\n\n";
#endif

    M_handler->addBC( name, flag, type, mode, base, comp );
}

// ===================================================
// Set Methods
// ===================================================
template< class Operator >
void BCInterface< Operator >::SetOperator( const boost::shared_ptr< Operator >& Oper )
{
    M_data.SetOperator( Oper );
}

template< class Operator >
void BCInterface< Operator >::SetHandler( const BCHandler_PtrType& handler )
{
    M_handler = handler;
}

template< class Operator >
void BCInterface< Operator >::SetHandlerParameters( const ID& bcNumber,
                                                    const BCHandler_Type::BCHints& hint )
{
    M_bcNumber = bcNumber;
    M_hint     = hint;

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::setHandlerParameters          M_bcNumber: " << M_bcNumber << "\n";
    Debug( 5020 ) << "                                               M_hint: " << M_hint << "\n";
#endif

}

// ===================================================
// Get Methods
// ===================================================
template< class Operator >
const typename BCInterface< Operator >::BCHandler_PtrType&
BCInterface< Operator >::GetHandler()
{
    return M_handler;
}

template< class Operator >
typename BCInterface< Operator >::Data_Type&
BCInterface< Operator >::GetDataContainer()
{
    return M_data;
}

// ===================================================
// Private Methods
// ===================================================
template< class Operator >
inline void
BCInterface< Operator >::BuildBase()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::BuildBase\n";
#endif

    //UInt position;

    switch ( M_data.GetBase().second )
    {
        case BCInterface_function:
        case BCInterface_functionFile:
        case BCInterface_OPERfunction:
        case BCInterface_OPERfunctionFile:

            AddBase( M_vectorFunction, M_data.GetBase().second );

            AddBCManager( M_vectorFunction.back()->GetBase() );

            break;

        case BCInterface_OPERFSI:

            AddBase( M_vectorFSI );

            AddBCManager( M_vectorFSI.back()->GetBase() );

            break;
    }
}

template< class Operator > template< class BCInterfaceBase >
inline void
BCInterface< Operator >::AddBase( std::vector< boost::shared_ptr< BCInterfaceBase > >& baseVector )
{
    boost::shared_ptr< BCInterfaceBase > Function( new BCInterfaceBase( M_data ) );
    baseVector.push_back( Function );
}

template< class Operator > template< class BCInterfaceBase >
inline void
BCInterface< Operator >::AddBase(       std::vector< boost::shared_ptr< BCInterfaceBase > >& baseVector,
                                  const BCBaseList_Type& Oper )
{
    boost::shared_ptr< BCInterfaceBase > Function( FactoryBCInterface_Function::instance().createObject( Oper ) );

    Function->SetData( M_data );

    baseVector.push_back( Function );
}

template< class Operator >
inline void
BCInterface< Operator >::AddBCManager( BCVectorInterface& base )
{
    if ( !M_handler.get() ) // If BCHandler has not been created yet, we do it now
        CreateHandler();

    switch ( M_data.GetMode() )
    {
        case Scalar:
        case Normal:
        case Tangential:
        case Directional:
        case Component:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5020 ) << "BCInterface::AddBCManager                              others" << "\n\n";
#endif

            std::cout << "ERROR: Scalar, Normal, Tangential, Directional, Component NOT AVAILABLE FOR FSI BC" << std::endl;
        break;

        case Full:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5020 ) << "BCInterface::AddBCManager                              Full" << "\n\n";
#endif

            M_handler->addBC( M_data.GetName(), M_data.GetFlag(), M_data.GetType(), M_data.GetMode(), base, M_data.GetComN() );

        break;
    }
}

template< class Operator > template< class BCBase >
inline void
BCInterface< Operator >::AddBCManager( BCBase& base )
{
    if ( !M_handler.get() ) // If BCHandler has not been created yet, we do it now
        CreateHandler();

    switch ( M_data.GetMode() )
    {
        case Scalar:
        case Normal:
        case Tangential:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5020 ) << "BCInterface::AddBCManager                              Scalar, Normal, Tangential" << "\n\n";
#endif

            M_handler->addBC( M_data.GetName(), M_data.GetFlag(), M_data.GetType(), M_data.GetMode(), base );

        break;

        case Directional:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5020 ) << "BCInterface::AddBCManager                              Directional" << "\n\n";
#endif
        {
            // Parameters for direction BC
            M_data.SetName( M_data.GetName() + "_direction" );
            M_data.SetBase( make_pair( "function", BCInterface_function ) );
            M_data.SetBaseString( M_data.GetDirection() );

            // Directional field
            AddBase( M_vectorFunction, M_data.GetBase().second );

            // Directional base
            boost::shared_ptr< BCFunctionDirectional > directionalBase( new BCFunctionDirectional( base.Function(), M_vectorFunction.back()->GetBase().Function() ) );
            M_vectorFunctionDirection.push_back( directionalBase );

            M_handler->addBC( M_data.GetName(), M_data.GetFlag(), M_data.GetType(), M_data.GetMode(), *M_vectorFunctionDirection.back() );
        }

        break;

        case Full:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5020 ) << "BCInterface::AddBCManager                              Full" << "\n\n";
#endif

            M_handler->addBC( M_data.GetName(), M_data.GetFlag(), M_data.GetType(), M_data.GetMode(), base, M_data.GetComN() );

        break;

        case Component:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5020 ) << "BCInterface::AddBCManager                              Component" << "\n\n";
#endif

            M_handler->addBC( M_data.GetName(), M_data.GetFlag(), M_data.GetType(), M_data.GetMode(), base, M_data.GetComV() );

        break;
    }
}

} // Namespace LifeV

#endif /* BCInterface_H */
