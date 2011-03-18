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
 *  @brief File containing the BCInterface class
 *
 *  @date 10-05-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface_H
#define BCInterface_H 1

#include <lifemc/lifesolver/BCInterfaceDefinitions.hpp>

#include <lifemc/lifesolver/BCInterfaceData.hpp>
#include <lifemc/lifesolver/BCInterfaceFactory.hpp>

namespace LifeV
{

//! BCInterface - LifeV Interface to load Boundary Conditions completely from a GetPot file
/*!
 *  @author Cristiano Malossi
 *
 *  This class allows to impose boundary conditions completely from a file.
 *
 *  <b>EXAMPLE - DATA FILE</b>
 *
 *  In the GetPot data file, BCInterface reads a new section: [boundary_conditions].
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
 *
 *  NOTE:
 *
 *  a) You can add manually more conditions by using setBC() after the call to buildHandler() function.
 *     In this case you have to manually set the TOTAL number of boundary conditions
 *     by using setHandlerParameters() function BEFORE building the handler.
 */
template< class BcHandler, class PhysicalSolverType >
class BCInterface
{
public:

    //! @name Type definitions
    //@{

    typedef BcHandler                                             bcHandler_Type;
    typedef boost::shared_ptr< bcHandler_Type >                   bcHandlerPtr_Type;

    typedef PhysicalSolverType                                    physicalSolver_Type;
    typedef boost::shared_ptr< physicalSolver_Type >              physicalSolverPtr_Type;

    typedef BCInterfaceFactory< physicalSolver_Type >             factory_Type;
    typedef typename factory_Type::bcFunctionPtr_Type             bcFunctionPtr_Type;

    typedef BCInterfaceData                                       data_Type;

    typedef std::vector< bcFunctionPtr_Type >                     vectorFunction_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface();

    //! Destructor
    virtual ~BCInterface() {}

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
        M_data.readBC( fileName, dataSection, name );
    }

    //! Insert the current boundary condition in the BChandler
    virtual void insertBC();

    //! Update the variables inside the physical solver
    void updatePhysicalSolverVariables();

    //@}


    //! @name Set Methods
    //@{

    //! Set a physical solver
    /*!
     * @param physicalSolver physical solver
     */
    virtual void setPhysicalSolver( const physicalSolverPtr_Type& physicalSolver );

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


protected:

    // Handler and parameters
    bcHandlerPtr_Type                        M_handler;

    // Data
    data_Type                                M_data;

    // Functions
    vectorFunction_Type                      M_vectorFunction;

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface( const BCInterface& interface1D );

    BCInterface& operator=( const BCInterface& interface1D );

    //@}
};

// ===================================================
// Constructors & Destructor
// ===================================================
template< class BcHandler, class PhysicalSolverType >
BCInterface< BcHandler, PhysicalSolverType >::BCInterface() :
        M_handler                 (),
        M_data                    (),
        M_vectorFunction          ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::BCInterface" << "\n";
#endif

}

// ===================================================
// Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >
void
BCInterface< BcHandler, PhysicalSolverType >::fillHandler( const std::string& fileName, const std::string& dataSection )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::fillHandler\n";
#endif

    GetPot dataFile( fileName );
    for ( UInt i( 0 ); i < dataFile.vector_variable_size( ( dataSection + "/boundary_conditions/list" ).c_str() ); ++i )
    {
        readBC( fileName, dataSection + "/boundary_conditions/",
                          dataFile( ( dataSection + "/boundary_conditions/list" ).c_str(), " ", i ) );

        this->insertBC();
    }
}

template< class BcHandler, class PhysicalSolverType >
inline void
BCInterface< BcHandler, PhysicalSolverType >::insertBC()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::insertBC\n";
#endif

    factory_Type factory;
    M_vectorFunction.push_back( factory.createFunction( M_data ) );
}

template< class BcHandler, class PhysicalSolverType >
void
BCInterface< BcHandler, PhysicalSolverType >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface::updatePhysicalSolverVariables\n";
#endif

    for ( UInt i( 0 ); i < M_vectorFunction.size(); ++i )
    {
        BCInterfaceFunctionSolver< physicalSolver_Type > *castedFunctionSolver =
            dynamic_cast < BCInterfaceFunctionSolver< physicalSolver_Type > * > ( &( *M_vectorFunction[i] ) );

        if ( castedFunctionSolver != 0 )
            castedFunctionSolver->updatePhysicalSolverVariables();
    }
}

// ===================================================
// Set Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >
void
BCInterface< BcHandler, PhysicalSolverType >::setPhysicalSolver( const physicalSolverPtr_Type& physicalSolver )
{
    //for ( typename vectorFunction_Type::const_iterator i = M_vectorFunction.begin() ; i < M_vectorFunction.end() ; ++i )
    for ( UInt i( 0 ); i < M_vectorFunction.size(); ++i )
    {
        BCInterfaceFunctionSolver< physicalSolver_Type > *castedFunctionSolver =
            dynamic_cast < BCInterfaceFunctionSolver< physicalSolver_Type > * > ( &( *M_vectorFunction[i] ) );

        if ( castedFunctionSolver != 0 )
            castedFunctionSolver->setPhysicalSolver( physicalSolver );
    }
}

} // Namespace LifeV

#endif /* BCInterface_H */
