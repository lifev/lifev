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

// BCInterface includes
#include <lifev/bc_interface/core/bc/BCInterface.hpp>
#include <lifev/bc_interface/1D/bc/BCInterfaceData1D.hpp>

// Template specializations
#include <lifev/bc_interface/1D/function/BCInterfaceFunctionParser1D.hpp>
#include <lifev/bc_interface/1D/function/BCInterfaceFunctionParserSolver1D.hpp>
#include <lifev/bc_interface/1D/function/BCInterfaceFunctionSolverDefined1D.hpp>
#include <lifev/bc_interface/1D/function/BCInterfaceFunctionUserDefined1D.hpp>

namespace LifeV
{

//! BCInterface1D - LifeV interface to load boundary conditions for 1D problems completely from a \c GetPot file
/*!
 *  @author Cristiano Malossi
 *
 *  This class allows to impose boundary conditions for a 1D problem completely from a file.
 *
 *  <b>EXAMPLE - DATA FILE</b> <BR>
 *  In the GetPot data file, \c BCInterface reads a new section: <CODE> [boundary_conditions] </CODE>.
 *
 *  Inside the new section there is a list of boundary conditions which correspond to other sub-section
 *  with the same name, for example: <CODE> list = 'InFlow OutFlow' </CODE>
 *
 *  Each boundary condition has a similar structure. The list of properties depends from the type of the
 *  boundary condition. For example:
 *
 *  <CODE>
 *  [InFlow]                             <BR>
 *  side                = left           <BR>
 *  quantity            = Q              <BR>
 *  line                = first          <BR>
 *  function            = 'sin(2*pi*t)'  <BR>
 *
 *  [OutFlow]                            <BR>
 *  side                = right          <BR>
 *  quantity            = W2             <BR>
 *  line                = first          <BR>
 *  functionSD          = Absorbing      <BR>
 *  </CODE>
 *
 *  where \c side, \c quantity, and \c line are the classical parameters for a 1D boundary condition.
 *  The string \c function represents the base module and can be replaced by other derived/alternative modules.
 *  The following functions are available (see the related classes for more information):
 *
 *  <ol>
 *      <li> \c function, which is implemented in \c BCInterfaceFunctionParser;
 *      <li> \c functionFile, which is implemented in \c BCInterfaceFunctionParserFile;
 *      <li> \c functionSolver, which is implemented in \c BCInterfaceFunctionParserSolver;
 *      <li> \c functionFileSolver, which is implemented in \c BCInterfaceFunctionParserFileSolver;
 *      <li> \c functionUD, which is implemented in \c BCInterfaceFunctionUserDefined;
 *      <li> \c functionSD, which is implemented in \c BCInterfaceFunctionSolverDefined.
 *  </ol>
 *
 *  All the parameters are case sensitive.
 *
 *  See \c BCInterface base class for more details.
 */
template< class BcHandler, class PhysicalSolverType >
class BCInterface1D : public virtual BCInterface< BcHandler, PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface< BcHandler, PhysicalSolverType >                bcInterface_Type;

    typedef typename bcInterface_Type::bcHandler_Type                   bcHandler_Type;
    typedef typename bcInterface_Type::bcHandlerPtr_Type                bcHandlerPtr_Type;

    typedef typename bcInterface_Type::physicalSolver_Type              physicalSolver_Type;
    typedef typename bcInterface_Type::physicalSolverPtr_Type           physicalSolverPtr_Type;

    typedef typename bcInterface_Type::factory_Type                     factory_Type;

    typedef typename bcInterface_Type::bcFunctionPtr_Type               bcFunctionPtr_Type;
    typedef typename bcInterface_Type::vectorFunction_Type              vectorFunction_Type;

    typedef typename bcInterface_Type::bcFunctionParserSolver_Type      bcFunctionParserSolver_Type;
    typedef typename bcInterface_Type::bcFunctionParserSolverPtr_Type   bcFunctionParserSolverPtr_Type;

    typedef typename bcInterface_Type::bcFunctionSolverDefinedPtr_Type  bcFunctionSolverDefinedPtr_Type;
    typedef typename bcInterface_Type::vectorFunctionSolverDefined_Type vectorFunctionSolverDefined_Type;

    typedef BCInterfaceData1D                                           data_Type;
    typedef boost::shared_ptr< data_Type >                              dataPtr_Type;

    typedef typename bcHandler_Type::solutionPtr_Type                   solutionPtr_Type;
    typedef typename bcHandler_Type::fluxPtr_Type                       fluxPtr_Type;
    typedef typename bcHandler_Type::sourcePtr_Type                     sourcePtr_Type;
    typedef typename bcHandler_Type::vectorPtrContainer_Type            vectorPtrContainer_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface1D() : bcInterface_Type(), M_data ( new data_Type() ) {}

    //! Destructor
    virtual ~BCInterface1D() {}

    //@}


    //! @name Methods
    //@{

    //! Read a specific boundary condition from a file and add it to the data container
    /*!
     * @param fileName Name of the data file
     * @param dataSection section in the data file
     * @param name name of the boundary condition
     */
    void readBC ( const std::string& fileName, const std::string& dataSection, const std::string& name )
    {
        M_data->readBC ( fileName, dataSection, name );
    }

    //! Insert the current boundary condition in the BChandler
    void insertBC();

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
    void setBC ( const OneDFSI::bcSide_Type& bcSide, const OneDFSI::bcLine_Type& bcLine, const OneDFSI::bcType_Type& bcType, const BCBaseType& base )
    {
        this->M_handler->setBC ( bcSide, bcLine, bcType, base );
    }

    //@}


    //! @name Set Methods
    //@{

    //! Set the solution for the members that need it
    /*!
     * @param flux flux
     * @param source source
     */
    void setFluxSource ( const fluxPtr_Type& flux, const sourcePtr_Type& source );

    //! Set the solution for the members that need it
    /*!
     * @param solution solution
     */
    void setSolution ( const solutionPtr_Type& solution );

    //@}


    //! @name Get Methods
    //@{

    //! Get the data container
    /*!
     * @return the data container
     */
    data_Type& dataContainer()
    {
        return *M_data;
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface1D ( const BCInterface1D& interface1D );

    BCInterface1D& operator= ( const BCInterface1D& interface1D );

    //@}


    //! @name Private Methods
    //@{

    template< class BCInterfaceBaseType >
    void createFunction ( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector );

    template< class BCBaseType >
    void addBcToHandler ( BCBaseType& base );

    //@}

    // Data
    dataPtr_Type                       M_data;
};

// ===================================================
// Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >
inline void
BCInterface1D< BcHandler, PhysicalSolverType >::insertBC()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5020 ) << "BCInterface1D::insertBC\n";
#endif

    // Definitions
    factory_Type factory;
    OneDFSIFunction base;

    // Define correct BCI type
    switch ( M_data->base().second )
    {
        case BCIFunctionParser:
        case BCIFunctionParserFile:
        case BCIFunctionParserSolver:
        case BCIFunctionParserFileSolver:
        case BCIFunctionUserDefined:
        {
            this->M_vectorFunction.push_back ( factory.createFunctionParser ( M_data ) );
            this->M_vectorFunction.back()->assignFunction ( base );

            break;
        }
        case BCIFunctionSolverDefined:
        {
            this->M_vectorFunctionSolverDefined.push_back ( factory.createFunctionSolverDefined ( M_data ) );
            this->M_vectorFunctionSolverDefined.back()->assignFunction ( base );

            break;
        }
        default:

            std::cout << " !!! Error: " << M_data->base().first << " is not valid in BCInterface1D !!!" << std::endl;
            return;
    }

    // Add base to BCHandler
    addBcToHandler ( base );
}

// ===================================================
// Set Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >
void
BCInterface1D< BcHandler, PhysicalSolverType >::setFluxSource ( const fluxPtr_Type& flux, const sourcePtr_Type& source )
{
    for ( typename vectorFunctionSolverDefined_Type::const_iterator i = this->M_vectorFunctionSolverDefined.begin() ; i < this->M_vectorFunctionSolverDefined.end() ; ++i )
    {
        ( *i )->setFluxSource ( flux, source );
    }

    this->M_handler->setFluxSource ( flux, source );
}

template< class BcHandler, class PhysicalSolverType >
void
BCInterface1D< BcHandler, PhysicalSolverType >::setSolution ( const solutionPtr_Type& solution )
{
    //for ( typename vectorFunction_Type::const_iterator i = M_vectorFunction.begin() ; i < M_vectorFunction.end() ; ++i )
    for ( UInt i ( 0 ); i < this->M_vectorFunction.size(); ++i )
    {
        bcFunctionParserSolverPtr_Type castedFunctionSolver = boost::dynamic_pointer_cast< bcFunctionParserSolver_Type > ( this->M_vectorFunction[i] );

        if ( castedFunctionSolver != 0 )
        {
            castedFunctionSolver->setSolution ( solution );
        }
    }

    for ( typename vectorFunctionSolverDefined_Type::const_iterator i = this->M_vectorFunctionSolverDefined.begin() ; i < this->M_vectorFunctionSolverDefined.end() ; ++i )
    {
        ( *i )->setSolution ( solution );
    }

    this->M_handler->setSolution ( solution );
}

// ===================================================
// Private Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType > template< class BCInterfaceBaseType >
inline void
BCInterface1D< BcHandler, PhysicalSolverType >::createFunction ( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector )
{
    boost::shared_ptr< BCInterfaceBaseType > function ( new BCInterfaceBaseType() );
    function->setData ( M_data );
    baseVector.push_back ( function );
}

template< class BcHandler, class PhysicalSolverType > template< class BCBaseType >
inline void
BCInterface1D< BcHandler, PhysicalSolverType >::addBcToHandler ( BCBaseType& base )
{
    if ( !this->M_handler.get() ) // If BCHandler has not been created yet, we do it now
    {
        this->createHandler();
    }

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5020 ) << "BCInterface1D::addBCManager" << "\n\n";
#endif

    this->M_handler->setBC ( M_data->side(), M_data->line(), M_data->quantity(), base );
}

} // Namespace LifeV

#endif /* BCInterface1D_H */
