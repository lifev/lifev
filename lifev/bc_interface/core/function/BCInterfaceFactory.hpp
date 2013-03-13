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
 *  @brief File containing the BCInterfaceFactory class
 *
 *  @date 18-03-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFactory_H
#define BCInterfaceFactory_H 1

#include <lifev/bc_interface/core/bc/BCInterfaceData.hpp>

#include <lifev/bc_interface/core/function/BCInterfaceFunctionParser.hpp>
#include <lifev/bc_interface/core/function/BCInterfaceFunctionParserFile.hpp>
#include <lifev/bc_interface/core/function/BCInterfaceFunctionParserSolver.hpp>
#include <lifev/bc_interface/core/function/BCInterfaceFunctionParserFileSolver.hpp>
#include <lifev/bc_interface/core/function/BCInterfaceFunctionUserDefined.hpp>
#include <lifev/bc_interface/core/function/BCInterfaceFunctionSolverDefined.hpp>

namespace LifeV
{

// Forward class declarations
template< typename BcHandlerType, typename PhysicalSolverType >
class BCInterfaceFunctionSolverDefined;

template< typename BcHandlerType, typename PhysicalSolverType >
inline BCInterfaceFunctionSolverDefined< BcHandlerType, PhysicalSolverType >* createBCInterfaceFunctionSolverDefined();



//! BCInterfaceFactory - Factory to create \c BCInterface functions
/*!
 *  @author Cristiano Malossi
 *
 *  This class allows to create boundary functions which can be used by any BCInterface implementation.
 *  The following functions are available (see the related classes for more information):
 *
 *  <ol>
 *      <li> \c function, which is implemented in \c BCInterfaceFunctionParser;
 *      <li> \c functionFile, which is implemented in \c BCInterfaceFunctionParserFile;
 *      <li> \c functionSolver, which is implemented in \c BCInterfaceFunctionParserSolver;
 *      <li> \c functionFileSolver, which is implemented in \c BCInterfaceFunctionParserFileSolver;
 *      <li> \c functionUD, which is implemented in \c BCInterfaceFunctionUserDefined;
 *      <li> \c functionSD, which is implemented in \c BCInterfaceFunctionSolverDefined;
 *  </ol>
 */

template< typename BcHandlerType, typename PhysicalSolverType >
class BCInterfaceFactory
{
public:

    //! @name Type definitions
    //@{

    typedef BcHandlerType                                                                                        bcHandler_Type;
    typedef PhysicalSolverType                                                                                   physicalSolver_Type;

    typedef BCInterfaceFunction< bcHandler_Type, physicalSolver_Type >                                           bcFunction_Type;
    typedef boost::shared_ptr< bcFunction_Type >                                                                 bcFunctionPtr_Type;
    typedef FactorySingleton< Factory< bcFunction_Type , baseList_Type > >                                       factoryFunction_Type;

    typedef BCInterfaceFunctionParserSolver< bcHandler_Type, physicalSolver_Type >                               bcFunctionParserSolver_Type;
    typedef boost::shared_ptr< bcFunctionParserSolver_Type >                                                     bcFunctionParserSolverPtr_Type;

    typedef BCInterfaceFunctionSolverDefined< bcHandler_Type, physicalSolver_Type >                              bcFunctionSolverDefined_Type;
    typedef boost::shared_ptr< bcFunctionSolverDefined_Type >                                                    bcFunctionSolverDefinedPtr_Type;
    typedef FactorySingleton< Factory< bcFunctionSolverDefined_Type, baseList_Type > >                           factoryFunctionSolverDefined_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceFactory();

    //! Destructor
    virtual ~BCInterfaceFactory() {}

    //@}


    //! @name Methods
    //@{

    //! Create a parser function
    /*!
     * @param data data container
     */
    template< typename DataType >
    bcFunctionPtr_Type createFunctionParser ( const DataType& data );

    //! Create a user defined function
    /*!
     * @param data data container
     */
    template< typename DataType >
    bcFunctionSolverDefinedPtr_Type createFunctionSolverDefined ( const DataType& data );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFactory ( const BCInterfaceFactory& bcInterfaceFactory );

    BCInterfaceFactory& operator= ( const BCInterfaceFactory& bcInterfaceFactory );

    //@}
};

// ===================================================
// Constructors & Destructor
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
BCInterfaceFactory< BcHandlerType, PhysicalSolverType >::BCInterfaceFactory()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5020 ) << "BCInterfaceFactory::BCInterfaceFactory" << "\n";
#endif

    //Factory registration
    factoryFunction_Type::instance().registerProduct (              BCIFunctionParser,           &createBCInterfaceFunctionParser< bcHandler_Type, physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct (              BCIFunctionParserFile,       &createBCInterfaceFunctionParserFile< bcHandler_Type, physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct (              BCIFunctionParserSolver,     &createBCInterfaceFunctionParserSolver< bcHandler_Type, physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct (              BCIFunctionParserFileSolver, &createBCInterfaceFunctionParserFileSolver< bcHandler_Type, physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct (              BCIFunctionUserDefined,      &createBCInterfaceFunctionUserDefined< bcHandler_Type, physicalSolver_Type > );
    factoryFunctionSolverDefined_Type::instance().registerProduct ( BCIFunctionSolverDefined,    &createBCInterfaceFunctionSolverDefined< bcHandler_Type, physicalSolver_Type > );
}

// ===================================================
// Methods
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType > template< typename DataType >
inline typename BCInterfaceFactory< BcHandlerType, PhysicalSolverType >::bcFunctionPtr_Type
BCInterfaceFactory< BcHandlerType, PhysicalSolverType >::createFunctionParser ( const DataType& data )
{
    bcFunctionPtr_Type function ( factoryFunction_Type::instance().createObject ( data.base().second, data.mapBase() ) );

    function->setData ( data );

    return function;
}

template< typename BcHandlerType, typename PhysicalSolverType > template< typename DataType >
inline typename BCInterfaceFactory< BcHandlerType, PhysicalSolverType >::bcFunctionSolverDefinedPtr_Type
BCInterfaceFactory< BcHandlerType, PhysicalSolverType >::createFunctionSolverDefined ( const DataType& data )
{
    bcFunctionSolverDefinedPtr_Type function ( factoryFunctionSolverDefined_Type::instance().createObject ( data.base().second, data.mapBase() ) );

    function->setData ( data );

    return function;
}

} // Namespace LifeV

#endif /* BCInterfaceFactory_H */
