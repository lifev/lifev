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
 *  @brief File containing the MultiScale Definitions
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

// Boost classes
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

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

#include <life/lifealg/EpetraMap.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <life/lifearray/EpetraMatrix.hpp>

namespace LifeV
{
namespace multiscale
{

// Enum objects
/*! @enum algorithmsTypes
 */
enum algorithms_Type
{
    Aitken,   /*!< Aitken method */
    Explicit, /*!< Explicit method */
    Newton    /*!< Newton method (with exact Jacobian Matrix) */
};

/*! @enum modelsTypes
 */
enum models_Type
{
    Fluid3D,                   /*!< Fluid (Oseen) 3D model */
    FSI3D,                     /*!< FSI 3D model */
    MultiScale,                /*!< MultiScale model */
    OneDimensional             /*!< 1D model */
};

/*! @enum couplingsTypes
 */
enum couplings_Type
{
    BoundaryCondition,         /*!< Boundary condition */
    FlowRateStress,            /*!< Flow Rate/stress coupling condition */
    Stress                     /*!< All stress coupling condition */
};

/*! @enum stressTypes
 */
enum stress_Type
{
    StaticPressure,             /*!< Use static pressure */
    TotalPressure,              /*!< Use total pressure (static + dynamic) */
    LagrangeMultiplier          /*!< Use the Lagrange multiplier */
};

enum errors_Type
{
    IterationsMaximumNumber, /*!< Maximum number of iterations reached */
    Tolerance,               /*!< Tolerance not satisfied */
    Residual,                /*!< External residual not satisfied */
    ModelType,               /*!< Model type not recognized */
    CouplingType             /*!< Coupling type not recognized */
};

// Folder of the problem
extern std::string multiscaleProblemFolder;

// Step of the problem ( > 0 when performing a restart )
extern UInt multiscaleProblemStep;

// Exit Flag
extern bool multiscaleExitFlag;

// Map objects
extern std::map< std::string, algorithms_Type > multiscaleAlgorithmsMap;
extern std::map< std::string, models_Type >     multiscaleModelsMap;
extern std::map< std::string, couplings_Type >  multiscaleCouplingsMap;
extern std::map< std::string, stress_Type >     multiscaleStressesMap;

// Forward class declarations
class MultiscaleAlgorithm;
class MultiscaleModel;
class MultiscaleCoupling;
class MultiscaleData;

// Type definitions
typedef EntityFlag                                                        BCFlag;

typedef Displayer::comm_PtrType                                           multiscaleCommPtr_Type;

typedef EpetraVector                                                      multiscaleVector_Type;
typedef boost::shared_ptr< multiscaleVector_Type >                        multiscaleVectorPtr_Type;

typedef EpetraMatrix< Real >                                              multiscaleMatrix_Type;
typedef boost::shared_ptr< multiscaleMatrix_Type >                        multiscaleMatrixPtr_Type;

typedef MultiscaleAlgorithm                                               multiscaleAlgorithm_Type;
typedef boost::shared_ptr< multiscaleAlgorithm_Type >                     multiscaleAlgorithmPtr_Type;
typedef singleton< factory< multiscaleAlgorithm_Type, algorithms_Type > > multiscaleAlgorithmFactory_Type;

typedef MultiscaleModel                                                   multiscaleModel_Type;
typedef boost::shared_ptr< multiscaleModel_Type >                         multiscaleModelPtr_Type;
typedef singleton< factory< multiscaleModel_Type, models_Type > >         multiscaleModelFactory_Type;

typedef MultiscaleCoupling                                                multiscaleCoupling_Type;
typedef boost::shared_ptr< multiscaleCoupling_Type >                      multiscaleCouplingPtr_Type;
typedef singleton< factory< multiscaleCoupling_Type, couplings_Type > >   multiscaleCouplingFactory_Type;

typedef std::vector< multiscaleModelPtr_Type >                            multiscaleModelsVector_Type;
typedef multiscaleModelsVector_Type::iterator                             multiscaleModelsVectorIterator_Type;
typedef multiscaleModelsVector_Type::const_iterator                       multiscaleModelsVectorConstIterator_Type;

typedef std::vector< multiscaleCouplingPtr_Type >                         multiscaleCouplingsVector_Type;
typedef multiscaleCouplingsVector_Type::iterator                          multiscaleCouplingsVectorIterator_Type;
typedef multiscaleCouplingsVector_Type::const_iterator                    multiscaleCouplingsVectorConstIterator_Type;

typedef MultiscaleData                                                    multiscaleData_Type;
typedef boost::shared_ptr< multiscaleData_Type >                          multiscaleDataPtr_Type;

// ===================================================
// MS Utility Methods
// ===================================================

//! Define the map of the MS objects
inline void
multiscaleMapsDefinition()
{
    multiscaleModelsMap["Fluid3D"]              = Fluid3D;
    multiscaleModelsMap["FSI3D"]                = FSI3D;
    multiscaleModelsMap["MultiScale"]           = MultiScale;
    multiscaleModelsMap["OneDimensional"]       = OneDimensional;

    multiscaleCouplingsMap["BoundaryCondition"] = BoundaryCondition;
    multiscaleCouplingsMap["FlowRateStress"]    = FlowRateStress;
    multiscaleCouplingsMap["Stress"]            = Stress;

    multiscaleAlgorithmsMap["Aitken"]           = Aitken;
    multiscaleAlgorithmsMap["Explicit"]         = Explicit;
    multiscaleAlgorithmsMap["Newton"]           = Newton;

    multiscaleStressesMap["StaticPressure"]     = StaticPressure;
    multiscaleStressesMap["TotalPressure"]      = TotalPressure;
    multiscaleStressesMap["LagrangeMultiplier"] = LagrangeMultiplier;
}

//! Perform a dynamic cast from a base class to a derived class
/*!
 * @param base - pointer to the base object
 * @return pointer to the derived object
 */
template < typename DerivedType, typename BasePtrType >
inline DerivedType*
multiscaleDynamicCast( BasePtrType& base )
{
    return dynamic_cast< DerivedType * > ( &( *base ) );
}

//! Display and error message
/*!
 * @param errorMessage - The message that should be displayed
 */
inline void
multiscaleErrorMessage( const std::stringstream& errorMessage )
{
    std::cout << "MS ERROR: " << errorMessage.str() << std::endl;
}

//! Create an error message
/*!
 * @param error - The error type
 * @param message - Additional information about the error
 */
inline void
multiscaleErrorCheck( const errors_Type& error, const std::string& message = "" )
{
    std::stringstream errorMessage;
    errorMessage << std::scientific << std::setprecision( 6 );

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

    case ModelType:

        errorMessage << "Model type incorrect!\n";

        break;

    case CouplingType:

        errorMessage << "Coupling type incorrect!\n";

        break;

    default:

        errorMessage << "No error message for this errorType!\n";
    }

    errorMessage << message << "\n";
    multiscaleErrorMessage( errorMessage );

    // Change ExitFlag
    multiscaleExitFlag = EXIT_FAILURE;
}

} // Namespace multiscale

} // Namespace LifeV

#endif /* MS_Definitions_H */
