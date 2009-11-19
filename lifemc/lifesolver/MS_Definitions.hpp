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
 *  @brief MultiScale Definitions
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 03-11-2009
 */

#ifndef MS_Definitions_H
#define MS_Definitions_H 1

//#define DEBUG 1;

// LifeV classes
#include <life/lifecore/life.hpp>
#include <life/lifecore/util_string.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/displayer.hpp>
#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>
#include <life/lifearray/tab.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifemesh/markers.hpp>

// LifeV Trilinos
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
    #include <Epetra_MpiComm.h>
#else
    #include <Epetra_SerialComm.h>
#endif

#include <life/lifealg/EpetraMap.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <life/lifearray/EpetraMatrix.hpp>

// Boost classes
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

// STL classes
#include <string>
#include <fstream>

namespace LifeV {

// Enum objects
/*! @enum algorithmsTypes
 */
enum algorithmsTypes
{
    Aitken, /*!< Aitken method */
    Newton  /*!< Newton method (with Jacobian Matrix)*/
};

/*! @enum modelsTypes
 */
enum modelsTypes
{
    MultiScale, /*!< MultiScale model */
    Fluid3D     /*!< Oseen fluid 3D model */
};

/*! @enum couplingsTypes
 */
enum couplingsTypes
{
    BoundaryCondition, /*!< Boundary condition */
    Stress,            /*!< All stess coupling condition */
    FluxStress         /*!< Flux/stress coupling condition */
};

/*! @enum stressTypes
 */
enum stressTypes
{
    StaticPressure,    /*!< Use static pressure */
    TotalPressure      /*!< Use tolal pressure (static + dynamic) */

};

enum errorsTypes
{
    MS_IterationsMaximumNumber, /*!< Maximum number of iterations reached */
    MS_Tolerance,               /*!< Tolerance not satisfied */
    MS_ModelType,               /*!< Model type not recognized */
    MS_CouplingType             /*!< Coupling type not recognized */
};

// Exit Flag
extern bool MS_ExitFlag;

// Map objects
extern std::map< std::string, algorithmsTypes > algorithmMap;
extern std::map< std::string, modelsTypes >     modelsMap;
extern std::map< std::string, couplingsTypes >  couplingsMap;
extern std::map< std::string, stressTypes >     stressMap;

// Forward class declarations
class MS_Algorithm;
class MS_PhysicalModel;
class MS_PhysicalCoupling;

typedef EntityFlag                                                BCFlag;

typedef EpetraVector                                              VectorType;
typedef boost::shared_ptr< VectorType >                           Vector_ptrType;

typedef EpetraMatrix< Real >                                      MatrixType;
typedef boost::shared_ptr< MatrixType >                           Matrix_ptrType;

typedef MS_Algorithm                                              AlgorithmType;
typedef boost::shared_ptr< AlgorithmType >                        Algorithm_ptrType;

typedef MS_PhysicalModel                                          ModelType;
typedef boost::shared_ptr< ModelType >                            Model_ptrType;

typedef MS_PhysicalCoupling                                       CouplingType;
typedef boost::shared_ptr< CouplingType >                         Coupling_ptrType;

typedef std::vector<Model_ptrType>                                ModelsVector_Type;
typedef ModelsVector_Type::iterator                               ModelsVector_Iterator;
typedef ModelsVector_Type::const_iterator                         ModelsVector_ConstIterator;

typedef std::vector<Coupling_ptrType>                             CouplingsVector_Type;
typedef CouplingsVector_Type::iterator                            CouplingsVector_Iterator;
typedef CouplingsVector_Type::const_iterator                      CouplingsVector_ConstIterator;

typedef singleton< factory< AlgorithmType, algorithmsTypes > >    FactoryAlgorithms;
typedef singleton< factory< ModelType, modelsTypes > >            FactoryModels;
typedef singleton< factory< CouplingType, couplingsTypes > >      FactoryCouplings;

// ===================================================
// MS Utility Methods
// ===================================================

//! Perform a dynamic cast from a base class to a derived class
/*!
 * @param base - pointer to the base object
 * @return pointer to the derived object
 */
template < typename DerivedType, typename Base_ptrType >
inline DerivedType*
MS_DynamicCast( Base_ptrType& base )
{
    return dynamic_cast< DerivedType * > ( &( *base ) );
}

//! Display and error message
/*!
 * @param errorMessage - The message that should be displayed
 */
inline void
MS_ErrorMessage( const std::stringstream& errorMessage )
{
    std::cout << "MS ERROR: " << errorMessage.str() << std::endl;
}

//! Create an error message
/*!
 * @param error - The error type
 * @param message - Additional information about the error
 */
inline void
MS_ErrorCheck( const errorsTypes& error, const std::string& message = "" )
{
    std::stringstream errorMessage;
    errorMessage << std::scientific << std::setprecision( 6 );

    switch ( error )
    {
        case MS_IterationsMaximumNumber:

            errorMessage << "Maximum number of iterations reached!\n";

            break;

        case MS_Tolerance:

            errorMessage << "Tolerance not satisfied!\n";

            break;

        case MS_ModelType:

            errorMessage << "Model type incorrect!\n";

            break;

        case MS_CouplingType:

            errorMessage << "Coupling type incorrect!\n";

            break;

        default:

            errorMessage << "No error message for this errorType!\n";
    }

    errorMessage << message << "\n";
    MS_ErrorMessage( errorMessage );

    // Change ExitFlag
    MS_ExitFlag = EXIT_FAILURE;
}

} // Namespace LifeV

#endif /* MS_Definitions_H */
