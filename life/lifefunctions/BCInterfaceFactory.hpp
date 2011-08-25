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

#include <life/lifesolver/BCInterfaceData.hpp>

#include <life/lifefunctions/BCInterfaceFunctionParser.hpp>
#include <life/lifefunctions/BCInterfaceFunctionParserFile.hpp>
#include <life/lifefunctions/BCInterfaceFunctionParserSolver.hpp>
#include <life/lifefunctions/BCInterfaceFunctionParserFileSolver.hpp>

#include <life/lifefunctions/BCInterfaceFunctionSolverDefined.hpp>

namespace LifeV
{

// Forward class declarations
template< class PhysicalSolverType >
class BCInterfaceFunctionSolverDefined;

template< typename PhysicalSolverType >
inline BCInterfaceFunctionSolverDefined< PhysicalSolverType >* createBCInterfaceFunctionSolverDefined();



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
 *      <li> \c functionSD, which is implemented in \c BCInterfaceFunctionSolverDefined;
 *  </ol>
 */

template< class PhysicalSolverType >
class BCInterfaceFactory
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                                                   physicalSolver_Type;

    typedef BCInterfaceFunctionParser< physicalSolver_Type >                                                     bcFunctionParser_Type;
    typedef boost::shared_ptr< bcFunctionParser_Type >                                                           bcFunctionParserPtr_Type;
    typedef FactorySingleton< Factory< bcFunctionParser_Type , baseList_Type > >                                 factoryFunctionParser_Type;

    typedef BCInterfaceFunctionSolverDefined< physicalSolver_Type >                                              bcFunctionSolverDefined_Type;
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
    bcFunctionParserPtr_Type createFunctionParser( const DataType& data );

    //! Create a user defined function
    /*!
     * @param data data container
     */
    template< typename DataType >
    bcFunctionSolverDefinedPtr_Type createFunctionSolverDefined( const DataType& data );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFactory( const BCInterfaceFactory& bcInterfaceFactory );

    BCInterfaceFactory& operator=( const BCInterfaceFactory& bcInterfaceFactory );

    //@}
};

// ===================================================
// Constructors & Destructor
// ===================================================
template< class PhysicalSolverType >
BCInterfaceFactory< PhysicalSolverType >::BCInterfaceFactory()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterfaceFactory::BCInterfaceFactory" << "\n";
#endif

    //Factory registration
    factoryFunctionParser_Type::instance().registerProduct( BCIFunctionParser,           &createBCInterfaceFunctionParser< physicalSolver_Type > );
    factoryFunctionParser_Type::instance().registerProduct( BCIFunctionParserFile,       &createBCInterfaceFunctionParserFile< physicalSolver_Type > );
    factoryFunctionParser_Type::instance().registerProduct( BCIFunctionParserSolver,     &createBCInterfaceFunctionParserSolver< physicalSolver_Type > );
    factoryFunctionParser_Type::instance().registerProduct( BCIFunctionParserFileSolver, &createBCInterfaceFunctionParserFileSolver< physicalSolver_Type > );

    factoryFunctionSolverDefined_Type::instance().registerProduct( BCIFunctionSolverDefined, &createBCInterfaceFunctionSolverDefined< physicalSolver_Type > );
}

// ===================================================
// Methods
// ===================================================
template< class PhysicalSolverType > template< typename DataType >
inline typename BCInterfaceFactory< PhysicalSolverType >::bcFunctionParserPtr_Type
BCInterfaceFactory< PhysicalSolverType >::createFunctionParser( const DataType& data )
{
    bcFunctionParserPtr_Type function( factoryFunctionParser_Type::instance().createObject( data.base().second, data.mapBase() ) );

    function->setData( data );

    return function;
}

template< class PhysicalSolverType > template< typename DataType >
inline typename BCInterfaceFactory< PhysicalSolverType >::bcFunctionSolverDefinedPtr_Type
BCInterfaceFactory< PhysicalSolverType >::createFunctionSolverDefined( const DataType& data )
{
    bcFunctionSolverDefinedPtr_Type function( factoryFunctionSolverDefined_Type::instance().createObject( data.base().second, data.mapBase() ) );

    function->setData( data );

    return function;
}

} // Namespace LifeV

#endif /* BCInterfaceFactory_H */
