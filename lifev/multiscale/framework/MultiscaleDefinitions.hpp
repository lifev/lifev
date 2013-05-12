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
 *  @brief File containing the Multiscale Definitions
 *
 *  @date 03-11-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MS_Definitions_H
#define MS_Definitions_H 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// STL classes
#include <string>
#include <fstream>
#include <sstream>

// Boost classes
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LifeV classes
#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/StringUtility.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <lifev/core/mesh/MarkerDefinitions.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/multiscale/framework/MultiscaleGlobalData.hpp>

namespace LifeV
{
namespace Multiscale
{

// Enum objects
/*! @enum algorithmsTypes
 */
enum algorithms_Type
{
    Aitken,                 /*!< Aitken method */
    Broyden,                /*!< Broyden method (start from exact Jacobian Matrix) */
    Explicit,               /*!< Explicit method */
    Newton                  /*!< Newton method (with exact Jacobian Matrix) */
};

/*! @enum modelsTypes
 */
enum models_Type
{
    Fluid3D,                /*!< Fluid (Oseen) 3D model */
    FSI1D,                  /*!< 1D model */
    FSI3D,                  /*!< FSI 3D model */
    Multiscale,             /*!< Multiscale model */
    Windkessel0D,           /*!< Windkessel0D model */
    ZeroDimensional         /*!< 0D model */
};

/*! @enum couplingsTypes
 */
enum couplings_Type
{
    BoundaryCondition,        /*!< Boundary condition */
    MeanNormalStress,         /*!< Mean normal stress coupling condition */
    MeanNormalStressArea,     /*!< Mean normal stress with area coupling condition */
    MeanNormalStressValve,    /*!< Mean normal stress coupling condition with simple valve*/
    MeanTotalNormalStress,    /*!< Mean total normal stress coupling condition */
    MeanTotalNormalStressArea /*!< Mean total normal stress with area coupling condition */
};

enum errors_Type
{
    IterationsMaximumNumber, /*!< Maximum number of iterations reached */
    Tolerance,               /*!< Tolerance not satisfied */
    Residual,                /*!< External residual not satisfied */
    Solution,                /*!< Solution check not satisfied */
    ModelType,               /*!< Model type not recognized */
    CouplingType,            /*!< Coupling type not recognized */
    ModelInterface           /*!< Model interface not available */
};

// Folder of the problem
extern UInt multiscaleCoresPerNode;

// Folder of the problem
extern std::string multiscaleProblemFolder;

// Prefix of the problem
extern std::string multiscaleProblemPrefix;

// Step of the problem ( > 0 when performing a restart )
extern UInt multiscaleProblemStep;

// Save each N time steps
extern UInt multiscaleSaveEachNTimeSteps;

// Exit Flag
extern bool multiscaleExitFlag;

// Map objects
extern std::map< std::string, algorithms_Type > multiscaleAlgorithmsMap;
extern std::map< std::string, models_Type >     multiscaleModelsMap;
extern std::map< std::string, couplings_Type >  multiscaleCouplingsMap;

// Forward class declarations
class MultiscaleAlgorithm;
class MultiscaleModel;
class MultiscaleCoupling;

// Type definitions
typedef flag_Type                                                                multiscaleID_Type;
typedef std::vector< multiscaleID_Type >                                         multiscaleIDContainer_Type;
typedef multiscaleIDContainer_Type::const_iterator                               multiscaleIDContainerConstIterator_Type;

typedef Displayer::commPtr_Type                                                  multiscaleCommPtr_Type;

typedef VectorEpetra                                                             multiscaleVector_Type;
typedef boost::shared_ptr< multiscaleVector_Type >                               multiscaleVectorPtr_Type;

typedef MatrixEpetra< Real >                                                     multiscaleMatrix_Type;
typedef boost::shared_ptr< multiscaleMatrix_Type >                               multiscaleMatrixPtr_Type;

typedef LinearSolver::parameterList_Type                                         multiscaleParameterList_Type;
typedef LinearSolver::parameterListPtr_Type                                      multiscaleParameterListPtr_Type;

typedef MultiscaleAlgorithm                                                      multiscaleAlgorithm_Type;
typedef boost::shared_ptr< multiscaleAlgorithm_Type >                            multiscaleAlgorithmPtr_Type;
typedef FactorySingleton< Factory< multiscaleAlgorithm_Type, algorithms_Type > > multiscaleAlgorithmFactory_Type;

typedef MultiscaleModel                                                          multiscaleModel_Type;
typedef boost::shared_ptr< multiscaleModel_Type >                                multiscaleModelPtr_Type;
typedef FactorySingleton< Factory< multiscaleModel_Type, models_Type > >         multiscaleModelFactory_Type;

typedef MultiscaleCoupling                                                       multiscaleCoupling_Type;
typedef boost::shared_ptr< multiscaleCoupling_Type >                             multiscaleCouplingPtr_Type;
typedef FactorySingleton< Factory< multiscaleCoupling_Type, couplings_Type > >   multiscaleCouplingFactory_Type;

typedef std::vector< multiscaleModelPtr_Type >                                   multiscaleModelsContainer_Type;
typedef multiscaleModelsContainer_Type::iterator                                 multiscaleModelsContainerIterator_Type;
typedef multiscaleModelsContainer_Type::const_iterator                           multiscaleModelsContainerConstIterator_Type;

typedef std::vector< multiscaleCouplingPtr_Type >                                multiscaleCouplingsContainer_Type;
typedef multiscaleCouplingsContainer_Type::iterator                              multiscaleCouplingsContainerIterator_Type;
typedef multiscaleCouplingsContainer_Type::const_iterator                        multiscaleCouplingsContainerConstIterator_Type;

typedef MultiscaleGlobalData                                                     multiscaleData_Type;
typedef boost::shared_ptr< multiscaleData_Type >                                 multiscaleDataPtr_Type;

// ===================================================
// Multiscale Utility Methods
// ===================================================

//! Define the map of the MS objects
inline void
multiscaleMapsDefinition()
{
#if defined(LIFEV_HAS_NAVIERSTOKES)
    multiscaleModelsMap["Fluid3D"]         = Fluid3D;
#endif
#if defined(LIFEV_HAS_ONEDFSI)
    multiscaleModelsMap["FSI1D"]           = FSI1D;
#endif
#if defined(LIFEV_HAS_FSI)
    multiscaleModelsMap["FSI3D"]           = FSI3D;
#endif
    multiscaleModelsMap["Multiscale"]      = Multiscale;
#if defined(LIFEV_HAS_ZERODIMENSIONAL)
    multiscaleModelsMap["Windkessel0D"]    = Windkessel0D;
    multiscaleModelsMap["ZeroDimensional"] = ZeroDimensional;
#endif

    multiscaleCouplingsMap["BoundaryCondition"]         = BoundaryCondition;
    multiscaleCouplingsMap["MeanNormalStress"]          = MeanNormalStress;
    multiscaleCouplingsMap["MeanNormalStressValve"]     = MeanNormalStressValve;
    multiscaleCouplingsMap["MeanTotalNormalStress"]     = MeanTotalNormalStress;
#if defined(LIFEV_HAS_ONEDFSI) && defined(LIFEV_HAS_FSI)
    multiscaleCouplingsMap["MeanNormalStressArea"]      = MeanNormalStressArea;
    multiscaleCouplingsMap["MeanTotalNormalStressArea"] = MeanTotalNormalStressArea;
#endif

    multiscaleAlgorithmsMap["Aitken"]   = Aitken;
    multiscaleAlgorithmsMap["Broyden"]  = Broyden;
    multiscaleAlgorithmsMap["Explicit"] = Explicit;
    multiscaleAlgorithmsMap["Newton"]   = Newton;
}

//! Perform a dynamic cast from a base class to a derived class
/*!
 * @param base - pointer to the base object
 * @return pointer to the derived object
 */
template < typename DerivedType, typename BasePtrType >
inline boost::shared_ptr< DerivedType >
multiscaleDynamicCast ( BasePtrType& base )
{
    return boost::dynamic_pointer_cast< DerivedType > ( base );
}

//! Display and error message
/*!
 * @param errorMessage - The message that should be displayed
 */
inline void
multiscaleErrorMessage ( const std::stringstream& errorMessage )
{
    std::cerr << std::setprecision ( 10 ) << std::scientific << "MS ERROR: " << errorMessage.str() << std::endl;
}

//! Create an error message
/*!
 * @param error - The error type
 * @param message - Additional information about the error
 */
inline void
multiscaleErrorCheck ( const errors_Type& error, const std::string& message = "", const UInt& isLeader = true )
{
    if ( isLeader )
    {
        std::stringstream errorMessage;

        switch ( error )
        {
            case IterationsMaximumNumber:

                errorMessage << "Maximum number of iterations reached!\n";

                break;

            case Tolerance:

                errorMessage << "Tolerance not satisfied!\n";

                break;

            case Residual:

                errorMessage << "External residual not satisfied!\n";

                break;

            case Solution:

                errorMessage << "Solution check not satisfied!\n";

                break;

            case ModelType:

                errorMessage << "Model type incorrect!\n";

                break;

            case CouplingType:

                errorMessage << "Coupling type incorrect!\n";

                break;

            case ModelInterface:

                errorMessage << "Model interface not available!\n";

                break;

            default:

                errorMessage << "No error message for this errorType!\n";

                break;
        }

        errorMessage << message << "\n";
        multiscaleErrorMessage ( errorMessage );
    }

    // Change ExitFlag
    multiscaleExitFlag = EXIT_FAILURE;
}

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MS_Definitions_H */
