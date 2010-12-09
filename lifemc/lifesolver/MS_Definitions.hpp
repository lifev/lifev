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

// Enum objects
/*! @enum algorithmsTypes
 */
enum algorithmsTypes
{
    Aitken,   /*!< Aitken method */
    Explicit, /*!< Explicit method */
    Newton    /*!< Newton method (with exact Jacobian Matrix) */
};

/*! @enum modelsTypes
 */
enum modelsTypes
{
    Fluid3D,                   /*!< Fluid (Oseen) 3D model */
    FSI3D,                     /*!< FSI 3D model */
    MultiScale,                /*!< MultiScale model */
    OneDimensional             /*!< 1D model */
};

/*! @enum couplingsTypes
 */
enum couplingsTypes
{
    BoundaryCondition,         /*!< Boundary condition */
    FluxStress,                /*!< Flux/stress coupling condition */
    Stress                     /*!< All stress coupling condition */
};

/*! @enum stressTypes
 */
enum stressTypes
{
    StaticPressure,             /*!< Use static pressure */
    TotalPressure,              /*!< Use total pressure (static + dynamic) */
    LagrangeMultiplier          /*!< Use the Lagrange multiplier */
};

enum errorsTypes
{
    MS_IterationsMaximumNumber, /*!< Maximum number of iterations reached */
    MS_Tolerance,               /*!< Tolerance not satisfied */
    MS_Residual,                /*!< External residual not satisfied */
    MS_ModelType,               /*!< Model type not recognized */
    MS_CouplingType             /*!< Coupling type not recognized */
};

// Folder of the problem
extern std::string MS_ProblemFolder;

// Step of the problem ( > 0 when performing a restart )
extern UInt MS_ProblemStep;

// Exit Flag
extern bool MS_ExitFlag;

// Map objects
extern std::map< std::string, algorithmsTypes > MS_algorithmsMap;
extern std::map< std::string, modelsTypes >     MS_modelsMap;
extern std::map< std::string, couplingsTypes >  MS_couplingsMap;
extern std::map< std::string, stressTypes >     MS_stressesMap;

// Forward class declarations
class MS_Algorithm;
class MS_PhysicalModel;
class MS_PhysicalCoupling;
class MS_PhysicalData;

// Type definitions
typedef EntityFlag                                                 BCFlag;

typedef Displayer::comm_PtrType                                    MS_Comm_PtrType;

typedef EpetraVector                                               MS_Vector_Type;
typedef boost::shared_ptr< MS_Vector_Type >                        MS_Vector_PtrType;

typedef EpetraMatrix< Real >                                       MS_Matrix_Type;
typedef boost::shared_ptr< MS_Matrix_Type >                        MS_Matrix_PtrType;

typedef MS_Algorithm                                               MS_Algorithm_Type;
typedef boost::shared_ptr< MS_Algorithm_Type >                     MS_Algorithm_PtrType;
typedef singleton< factory< MS_Algorithm_Type, algorithmsTypes > > MS_Algorithm_Factory;

typedef MS_PhysicalModel                                           MS_Model_Type;
typedef boost::shared_ptr< MS_Model_Type >                         MS_Model_PtrType;
typedef singleton< factory< MS_Model_Type, modelsTypes > >         MS_Model_Factory;

typedef MS_PhysicalCoupling                                        MS_Coupling_Type;
typedef boost::shared_ptr< MS_Coupling_Type >                      MS_Coupling_PtrType;
typedef singleton< factory< MS_Coupling_Type, couplingsTypes > >   MS_Coupling_Factory;

typedef std::vector< MS_Model_PtrType >                            MS_ModelsVector_Type;
typedef MS_ModelsVector_Type::iterator                             MS_ModelsVector_Iterator;
typedef MS_ModelsVector_Type::const_iterator                       MS_ModelsVector_ConstIterator;

typedef std::vector< MS_Coupling_PtrType >                         MS_CouplingsVector_Type;
typedef MS_CouplingsVector_Type::iterator                          MS_CouplingsVector_Iterator;
typedef MS_CouplingsVector_Type::const_iterator                    MS_CouplingsVector_ConstIterator;

typedef MS_PhysicalData                                            MS_GlobalDataContainer_Type;
typedef boost::shared_ptr< MS_GlobalDataContainer_Type >           MS_GlobalDataContainer_PtrType;

// ===================================================
// MS Utility Methods
// ===================================================

//! Define the map of the MS objects
inline void
MS_MapsDefinition()
{
    MS_modelsMap["Fluid3D"]              = Fluid3D;
    MS_modelsMap["FSI3D"]                = FSI3D;
    MS_modelsMap["MultiScale"]           = MultiScale;
    MS_modelsMap["OneDimensional"]       = OneDimensional;

    MS_couplingsMap["BoundaryCondition"] = BoundaryCondition;
    MS_couplingsMap["FluxStress"]        = FluxStress;
    MS_couplingsMap["Stress"]            = Stress;

    MS_algorithmsMap["Aitken"]           = Aitken;
    MS_algorithmsMap["Explicit"]         = Explicit;
    MS_algorithmsMap["Newton"]           = Newton;

    MS_stressesMap["StaticPressure"]     = StaticPressure;
    MS_stressesMap["TotalPressure"]      = TotalPressure;
    MS_stressesMap["LagrangeMultiplier"] = LagrangeMultiplier;
}

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

    case MS_Residual:

        errorMessage << "External residual not satisfied!\n";

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
