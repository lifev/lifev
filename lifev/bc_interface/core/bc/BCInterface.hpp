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
 *  @brief File containing the BCInterface main class
 *
 *  @date 10-05-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface_H
#define BCInterface_H 1

// BCInterface includes
#include <lifev/bc_interface/core/bc/BCInterfaceData.hpp>
#include <lifev/bc_interface/core/function/BCInterfaceFactory.hpp>

namespace LifeV
{

//! BCInterface - LifeV interface to load boundary conditions completely from a \c GetPot file
/*!
 *  @author Cristiano Malossi
 *
 *  This class allows to impose boundary conditions completely from a \c GetPot file. You can derive a
 *  specific implementation for the BCHandler of your problem. For 0D (\c BCInterface0D), 1D (\c BCInterface1D) and 3D (\c BCInterface3D) problems the derived classes
 *  are already available.
 *
 *  <b>EXAMPLE - DATA FILE</b> <BR>
 *  In the GetPot data file, \c BCInterface reads a new section: <CODE> [boundary_conditions] </CODE>.
 *
 *  Inside the new section there is a list of boundary conditions which correspond to other sub-section
 *  with the same name, for example: <CODE> list = 'InFlow OutFlow' </CODE>
 *
 *  Each boundary condition has a similar structure. The list of properties depends from the type of the
 *  boundary condition (if it is for 1D or 3D solver for example). For example:
 *
 *  <CODE>
 *  [InFlow]                             <BR>
 *  ...                                  <BR>
 *  ...                                  <BR>
 *  function   = '3*0.03*(1/4-(x^2+y^2)' <BR>
 *
 *  [OutFlow]        <BR>
 *  ...              <BR>
 *  ...              <BR>
 *  function   = '0' <BR>
 *  </CODE>
 *
 *  The string \c function represents the base module and can be replaced by other derived/alternative modules.
 *  The following functions are available (see the related classes for more information):
 *
 *  <ol>
 *      <li> \c function, which is implemented in \c BCInterfaceFunctionParser;
 *      <li> \c functionFile, which is implemented in \c BCInterfaceFunctionParserFile;
 *      <li> \c functionSolver, which is implemented in \c BCInterfaceFunctionParserSolver;
 *      <li> \c functionFileSolver, which is implemented in \c BCInterfaceFunctionParserFileSolver;
 *  </ol>
 *
 *  All the parameters are case sensitive.
 *
 *  <b>EXAMPLE - HOW TO USE</b> <BR>
 *  Here there is a short guide on how to create and use a BCInterface object.
 *
 *  <ol>
 *      <li> First of all, you have to define a BCInterface class:
 *
 *      <CODE>
 *      BCInterface bcInterface;
 *      </CODE>
 *
 *      <li> You can create an empty handler by calling:
 *
 *      <CODE>
 *      bcInterface.createHandler()
 *      </CODE>
 *
 *      or you can set it from outside
 *
 *      <CODE>
 *      std::shared_ptr< bcHandler_Type > bcHandler( new bcHandler_Type() ); <BR>
 *      bcInterface.setHandler( bcHandler );
 *      </CODE>
 *
 *      <li> Then you can add all the file boundary condition by calling
 *
 *      <CODE>
 *      bcInterface.fillHandler( "fileName.dat", "section" );
 *      </CODE>
 *
 *      Or you can add one specific boundary conditions by calling
 *
 *      <CODE>
 *      bcInterface.readBC( "fileName.dat", "section", "bcSection" ); <BR>
 *      bcInterface.insertBC();
 *      </CODE>
 *
 *      Note that between readBC and insertBC you can manipulate the BC parameters by accessing
 *      the data container:
 *
 *      <CODE>
 *      bcInterface.dataContainer();
 *      </CODE>
 *
 *      In addition, you can also add BC directly from the code by accessing the bcHandler
 *
 *      <CODE>
 *      M_bc.handler()->addBC( ... );
 *      </CODE>
 *
 *      <li> If you are using functions that use solver variables first you have to pass the solver
 *
 *      <CODE>
 *      M_bc.setPhysicalSolver( physicalSolverPtr );
 *      </CODE>
 *
 *      Then be sure to update the variable at each time step before using the BCHandler:
 *
 *      <CODE>
 *      M_bc.updatePhysicalSolverVariables();
 *      </CODE>
 *
 *      <li> Finally, to get the handler you can use:
 *
 *      <CODE>
 *      M_bc.handler();
 *      </CODE>
 *
 *  </ol>
 */
template< class BcHandler, class PhysicalSolverType >
class BCInterface
{
public:

    //! @name Type definitions
    //@{

    typedef BcHandler                                                 bcHandler_Type;
    typedef std::shared_ptr< bcHandler_Type >                       bcHandlerPtr_Type;

    typedef PhysicalSolverType                                        physicalSolver_Type;
    typedef std::shared_ptr< physicalSolver_Type >                  physicalSolverPtr_Type;

    typedef BCInterfaceFactory< bcHandler_Type, physicalSolver_Type > factory_Type;

    typedef typename factory_Type::bcFunctionPtr_Type                 bcFunctionPtr_Type;
    typedef std::vector< bcFunctionPtr_Type >                         vectorFunction_Type;

    typedef typename factory_Type::bcFunctionParserSolver_Type        bcFunctionParserSolver_Type;
    typedef typename factory_Type::bcFunctionParserSolverPtr_Type     bcFunctionParserSolverPtr_Type;

    typedef typename factory_Type::bcFunctionSolverDefinedPtr_Type    bcFunctionSolverDefinedPtr_Type;
    typedef std::vector< bcFunctionSolverDefinedPtr_Type >            vectorFunctionSolverDefined_Type;

    typedef BCInterfaceData                                           data_Type;
    typedef std::shared_ptr< data_Type >                            dataPtr_Type;

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
    void createHandler()
    {
        M_handler.reset ( new bcHandler_Type() );
    }

    //! Fill the bcHandler with the BC provided in the file.
    /*!
     * @param fileName Name of the data file
     * @param dataSection Subsection inside [boundary_conditions]
     */
    void fillHandler ( const std::string& fileName, const std::string& dataSection );

    //! Read a specific boundary condition from a file and add it to the data container
    /*!
     * @param fileName Name of the data file
     * @param dataSection section in the data file
     * @param name name of the boundary condition
     */
    virtual void readBC ( const std::string& fileName, const std::string& dataSection, const std::string& name ) = 0;

    //! Insert the current boundary condition in the BChandler
    virtual void insertBC() = 0;

    //! Update the variables inside the physical solver
    virtual void updatePhysicalSolverVariables();

    //@}


    //! @name Set Methods
    //@{

    //! Set a physical solver
    /*!
     * @param physicalSolver physical solver
     */
    virtual void setPhysicalSolver ( const physicalSolverPtr_Type& physicalSolver );

    //! Set an Handler
    /*!
     * @param handler BCHandler
     */
    void setHandler ( const bcHandlerPtr_Type& handler )
    {
        M_handler = handler;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the shared_ptr to the BCHandler
    /*!
     * @return the pointer to the BCHandler
     */
    bcHandlerPtr_Type& handler()
    {
        return M_handler;
    }

    //! Get the data container
    /*!
     * @return the data container
     */
    virtual data_Type& dataContainer() = 0;

    //@}


protected:

    // Handler and parameters
    bcHandlerPtr_Type                        M_handler;

    // Parser functions
    vectorFunction_Type                      M_vectorFunction;

    // User defined functions
    vectorFunctionSolverDefined_Type         M_vectorFunctionSolverDefined;

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface ( const BCInterface& interface1D );

    BCInterface& operator= ( const BCInterface& interface1D );

    //@}
};

// ===================================================
// Constructors & Destructor
// ===================================================
template< class BcHandler, class PhysicalSolverType >
BCInterface< BcHandler, PhysicalSolverType >::BCInterface() :
    M_handler                     (),
    M_vectorFunction              (),
    M_vectorFunctionSolverDefined ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5020 ) << "BCInterface::BCInterface" << "\n";
#endif

}

// ===================================================
// Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >
void
BCInterface< BcHandler, PhysicalSolverType >::fillHandler ( const std::string& fileName, const std::string& dataSection )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5020 ) << "BCInterface::fillHandler\n";
#endif

    GetPot dataFile ( fileName );
    for ( UInt i ( 0 ); i < dataFile.vector_variable_size ( ( dataSection + "/boundary_conditions/list" ).c_str() ); ++i )
    {
        readBC ( fileName, dataSection + "/boundary_conditions/",
                 dataFile ( ( dataSection + "/boundary_conditions/list" ).c_str(), " ", i ) );

        this->insertBC();
    }
}

template< class BcHandler, class PhysicalSolverType >
void
BCInterface< BcHandler, PhysicalSolverType >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5020 ) << "BCInterface::updatePhysicalSolverVariables\n";
#endif

    for ( UInt i ( 0 ); i < M_vectorFunction.size(); ++i )
    {
        bcFunctionParserSolverPtr_Type castedFunctionSolver = std::dynamic_pointer_cast< bcFunctionParserSolver_Type > ( M_vectorFunction[i] );

        if ( castedFunctionSolver != 0 )
        {
            castedFunctionSolver->updatePhysicalSolverVariables();
        }
    }

    for ( typename vectorFunctionSolverDefined_Type::const_iterator i = M_vectorFunctionSolverDefined.begin() ; i < M_vectorFunctionSolverDefined.end() ; ++i )
    {
        ( *i )->updatePhysicalSolverVariables();
    }
}

// ===================================================
// Set Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >
void
BCInterface< BcHandler, PhysicalSolverType >::setPhysicalSolver ( const physicalSolverPtr_Type& physicalSolver )
{
    //for ( typename vectorFunction_Type::const_iterator i = M_vectorFunction.begin() ; i < M_vectorFunction.end() ; ++i )
    for ( UInt i ( 0 ); i < M_vectorFunction.size(); ++i )
    {
        bcFunctionParserSolverPtr_Type castedFunctionSolver = std::dynamic_pointer_cast< bcFunctionParserSolver_Type > ( M_vectorFunction[i] );

        if ( castedFunctionSolver != 0 )
        {
            castedFunctionSolver->setPhysicalSolver ( physicalSolver );
        }
    }

    for ( typename vectorFunctionSolverDefined_Type::const_iterator i = M_vectorFunctionSolverDefined.begin() ; i < M_vectorFunctionSolverDefined.end() ; ++i )
    {
        ( *i )->setPhysicalSolver ( physicalSolver );
    }
}

} // Namespace LifeV

#endif /* BCInterface_H */
