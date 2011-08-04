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

#include <life/lifesolver/BCInterfaceDefinitions.hpp>

#include <life/lifesolver/BCInterfaceData.hpp>
#include <life/lifesolver/BCInterfaceFunction.hpp>
#include <life/lifesolver/BCInterfaceFunctionFile.hpp>
#include <life/lifesolver/BCInterfaceFunctionSolver.hpp>
#include <life/lifesolver/BCInterfaceFunctionFileSolver.hpp>

namespace LifeV
{

//! BCInterfaceFactory - Factory to create \c BCInterface functions
/*!
 *  @author Cristiano Malossi
 *
 *  This class allows to create boundary functions which can be used by any BCInterface implementation.
 *  The following functions are available (see the related classes for more information):
 *
 *  <ol>
 *      <li> \c function, which is implemented in \c BCInterfaceFunction;
 *      <li> \c functionFile, which is implemented in \c BCInterfaceFunctionFile;
 *      <li> \c functionSolver, which is implemented in \c BCInterfaceFunctionSolver;
 *      <li> \c functionFileSolver, which is implemented in \c BCInterfaceFunctionFileSolver;
 *  </ol>
 */

template< class PhysicalSolverType >
class BCInterfaceFactory
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                                        physicalSolver_Type;
    typedef BCInterfaceData                                                                           data_Type;

    typedef FactorySingleton< Factory< BCInterfaceFunction< physicalSolver_Type > , baseList_Type > > factoryFunction_Type;

    typedef BCInterfaceFunction< physicalSolver_Type >                                                bcFunction_Type;
    typedef boost::shared_ptr< bcFunction_Type >                                                      bcFunctionPtr_Type;

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

    //! Create a function
    /*!
     * @param data data container
     */
    bcFunctionPtr_Type createFunction( const data_Type& data );

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
    factoryFunction_Type::instance().registerProduct( BCIFunction,           &createBCInterfaceFunction< physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct( BCIFunctionFile,       &createBCInterfaceFunctionFile< physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct( BCIFunctionSolver,     &createBCInterfaceFunctionSolver< physicalSolver_Type > );
    factoryFunction_Type::instance().registerProduct( BCIFunctionFileSolver, &createBCInterfaceFunctionFileSolver< physicalSolver_Type > );
}

// ===================================================
// Methods
// ===================================================
template< class PhysicalSolverType >
inline typename BCInterfaceFactory< PhysicalSolverType >::bcFunctionPtr_Type
BCInterfaceFactory< PhysicalSolverType >::createFunction( const data_Type& data )
{
    bcFunctionPtr_Type function( factoryFunction_Type::instance().createObject( data.base().second, data.mapBase() ) );

    function->setData( data );

    return function;
}

} // Namespace LifeV

#endif /* BCInterfaceFactory_H */
