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
 *  @brief BCInterface1D
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 10-05-2010
 */

#ifndef BCInterface1D_H
#define BCInterface1D_H 1

#include <lifemc/lifesolver/BCInterface1D_Definitions.hpp>

#include <lifemc/lifesolver/BCInterface1D_Data.hpp>
#include <lifemc/lifesolver/BCInterface1D_Function.hpp>
#include <lifemc/lifesolver/BCInterface1D_FunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface1D_OperatorFunction.hpp>
#include <lifemc/lifesolver/BCInterface1D_OperatorFunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface1D_DefaultFunctions.hpp>

namespace LifeV {

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
 *  a) You can add manually more conditions by using setBC() after the call to buildHandler() function.
 *     In this case you have to manually set the TOTAL number of boundary conditions
 *     by using setHandlerParameters() function BEFORE building the handler.
 */
template< class Operator = OneDimensionalModel_Solver >
class BCInterface1D
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface1D_BaseList                                                         BCBaseList_Type;

    typedef singleton< factory< BCInterface1D_Function< Operator > , BCBaseList_Type > >   FactoryBCInterface_Function;

    typedef OneDimensionalModel_BCHandler                                                  BCHandler_Type;
    typedef boost::shared_ptr< BCHandler_Type >                                            BCHandler_PtrType;

    typedef BCHandler_Type::Solution_PtrType                                               Solution_PtrType;

    typedef BCInterface1D_Data< Operator >                                                 Data_Type;

    typedef std::vector< boost::shared_ptr< BCInterface1D_Function< Operator > > >         VectorFunction_Type;
    typedef std::vector< boost::shared_ptr< BCInterface1D_DefaultFunctions< Operator > > > VectorDefaultFunction_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface1D();

    //! Copy constructor
    /*!
     * @param interface	- BCInterface1D
     */
    BCInterface1D( const BCInterface1D& interface );

    //! Destructor
    ~BCInterface1D() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param interface BCInterface1D
     * @return reference to a copy of the class
     */
    BCInterface1D& operator=( const BCInterface1D& interface );

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
     * @param side side of the condition
     * @param line line of the condition
     * @param type type of the condition
     * @param base base of the condition
     */
    template< class BCBase >
    void setBC( const OneD_BCSide& bcSide,
                const OneD_BCLine& bcLine,
                const OneD_BC&     bcType,
                const BCBase&      base );

    //@}


    //! @name Set Methods
    //@{

    //! Set an operator
    /*!
     * @param Oper operator
     */
    void SetOperator( const boost::shared_ptr< Operator >& Oper );

    //! Set the solution for the members that need it
    /*!
     * @param solution solution
     */
    void SetSolution( const Solution_PtrType solution );

    //! Set an Handler
    /*!
     * @param handler BCHandler
     */
    void SetHandler( const BCHandler_PtrType& handler );

    //@}


    //! @name Get Methods
    //@{

    //! Get the shared_ptr to the BCHandler
    /*!
     * @return the pointer to the BCHandler
     */
    const BCHandler_PtrType& GetHandler() const;

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

    template< class BCBase >
    inline void AddBCManager( BCBase& base );

    //@}

    // Handler and parameters
    BCHandler_PtrType                        M_handler;

    // Data
    Data_Type                                M_data;

    // Functions
    VectorFunction_Type                      M_vectorFunction;

    // Default Functions
    VectorDefaultFunction_Type               M_vectorDefaultFunction1D;
};

// ===================================================
// Constructors & Destructor
// ===================================================
template< class Operator >
BCInterface1D< Operator >::BCInterface1D( ) :
    M_handler                 (),
    M_data                    (),
    M_vectorFunction          (),
    M_vectorDefaultFunction1D ()
{

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface1D::BCInterface1D------------------------------" << "\n";
#endif

    //Factory registration
    FactoryBCInterface_Function::instance().registerProduct( BCInterface1D_function,         &BCInterface1D_CreateFunction< Operator > );
    FactoryBCInterface_Function::instance().registerProduct( BCInterface1D_functionFile,     &BCInterface1D_CreateFunctionFile< Operator > );
    FactoryBCInterface_Function::instance().registerProduct( BCInterface1D_OPERfunction,     &BCInterface1D_CreateOperatorFunction< Operator > );
    FactoryBCInterface_Function::instance().registerProduct( BCInterface1D_OPERfunctionFile, &BCInterface1D_CreateOperatorFunctionFile< Operator > );
}

template< class Operator >
BCInterface1D< Operator >::BCInterface1D( const BCInterface1D& interface ) :
    M_handler                 ( interface.M_handler ),
    M_data                    ( interface.M_data ),
    M_vectorFunction          ( interface.M_vectorFunction ),
    M_vectorDefaultFunction1D ( interface.M_vectorDefaultFunction1D )
{
}

// ===================================================
// Operators
// ===================================================
template< class Operator >
BCInterface1D< Operator >&
BCInterface1D< Operator >::operator=( const BCInterface1D& interface )
{
    if ( this != &interface )
    {
        M_handler                 = interface.M_handler;
        M_data                    = interface.M_data;
        M_vectorFunction          = interface.M_vectorFunction;
        M_vectorDefaultFunction1D = interface.M_vectorDefaultFunction1D;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
template< class Operator >
void
BCInterface1D< Operator >::CreateHandler()
{

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface1D::CreateHandler\n";
#endif

    M_handler.reset( new BCHandler_Type() );
}

template< class Operator >
void
BCInterface1D< Operator >::FillHandler( const std::string& FileName,
                                        const std::string& dataSection )
{

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface1D::buildHandler\n";
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
BCInterface1D< Operator >::ReadBC( const std::string& FileName,
                                 const std::string& dataSection,
                                 const BCName&      name )
{

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface1D::ReadBC\n";
#endif

    M_data.ReadBC( FileName, dataSection, name );
}

template< class Operator >
void
BCInterface1D< Operator >::InsertBC()
{

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface1D::InsertBC\n";
#endif

    BuildBase();
}

template< class Operator >
void
BCInterface1D< Operator >::UpdateOperatorVariables()
{

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface1D::UpdateOperatorVariables\n";
#endif

    for ( UInt i( 0 ); i < M_vectorFunction.size(); ++i )
    {
        BCInterface1D_OperatorFunction< Operator > *Oper =
                dynamic_cast< BCInterface1D_OperatorFunction< Operator > * > ( &( *M_vectorFunction[i] ) );

        if ( Oper != 0 )
            Oper->UpdateOperatorVariables();
    }
}

template< class Operator > template< class BCBase >
void
BCInterface1D< Operator >::setBC( const OneD_BCSide& bcSide,
                                  const OneD_BCLine& bcLine,
                                  const OneD_BC&     bcType,
                                  const BCBase&      base )
{

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface1D::setBC" << "\n\n";
#endif

    M_handler->setBC( bcSide, bcLine, bcType, base );
}

// ===================================================
// Set Methods
// ===================================================
template< class Operator >
void BCInterface1D< Operator >::SetOperator( const boost::shared_ptr< Operator >& Oper )
{
    M_data.SetOperator( Oper );
}

template< class Operator >
void BCInterface1D< Operator >::SetSolution( const Solution_PtrType solution )
{
    //for ( typename VectorFunction_Type::const_iterator i = M_vectorFunction.begin() ; i < M_vectorFunction.end() ; ++i )
    for ( UInt i( 0 ); i < M_vectorFunction.size(); ++i )
    {
        BCInterface1D_OperatorFunction< Operator > *Oper =
                dynamic_cast < BCInterface1D_OperatorFunction< Operator >* > ( &( *M_vectorFunction[i] ) );

        if ( Oper != 0 )
            Oper->SetSolution( solution );
    }

    for ( typename VectorDefaultFunction_Type::const_iterator i = M_vectorDefaultFunction1D.begin() ; i < M_vectorDefaultFunction1D.end() ; ++i )
        ( *i )->SetSolution( solution );

    M_handler->setSolution( solution );
}

template< class Operator >
void BCInterface1D< Operator >::SetHandler( const BCHandler_PtrType& handler )
{
    M_handler = handler;
}

// ===================================================
// Get Methods
// ===================================================
template< class Operator >
const typename BCInterface1D< Operator >::BCHandler_PtrType&
BCInterface1D< Operator >::GetHandler() const
{
    return M_handler;
}

template< class Operator >
typename BCInterface1D< Operator >::Data_Type&
BCInterface1D< Operator >::GetDataContainer()
{
    return M_data;
}

// ===================================================
// Private Methods
// ===================================================
template< class Operator >
inline void
BCInterface1D< Operator >::BuildBase()
{

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface1D::BuildBase\n";
#endif

    //UInt position;

    switch ( M_data.GetBase().second )
    {
        case BCInterface1D_function:
        case BCInterface1D_functionFile:
        case BCInterface1D_OPERfunction:
        case BCInterface1D_OPERfunctionFile:

            AddBase( M_vectorFunction, M_data.GetBase().second );

            AddBCManager( M_vectorFunction.back()->GetBase() );

            break;

        case BCInterface1D_Default:

            AddBase( M_vectorDefaultFunction1D );

            AddBCManager( M_vectorDefaultFunction1D.back()->GetBase() );

            break;
    }
}

template< class Operator > template< class BCInterfaceBase >
inline void
BCInterface1D< Operator >::AddBase( std::vector< boost::shared_ptr< BCInterfaceBase > >& baseVector )
{
    boost::shared_ptr< BCInterfaceBase > Function( new BCInterfaceBase( M_data ) );
    baseVector.push_back( Function );
}

template< class Operator > template< class BCInterfaceBase >
inline void
BCInterface1D< Operator >::AddBase(       std::vector< boost::shared_ptr< BCInterfaceBase > >& baseVector,
                                    const BCBaseList_Type& Oper )
{
    boost::shared_ptr< BCInterfaceBase > Function( FactoryBCInterface_Function::instance().createObject( Oper ) );

    Function->SetData( M_data );

    baseVector.push_back( Function );
}

template< class Operator > template< class BCBase >
inline void
BCInterface1D< Operator >::AddBCManager( BCBase& base )
{
    if ( !M_handler.get() ) // If BCHandler has not been created yet, we do it now
        CreateHandler();

    M_handler->setBC( M_data.GetSide(), M_data.GetLine(), M_data.GetQuantity(), base );
}

} // Namespace LifeV

#endif /* BCInterface1D_H */
