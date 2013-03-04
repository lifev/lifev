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
 *  @brief One Dimensional Model Global Definitions
 *
 *  @version 1.0
 *  @date 15-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDFSIDefinitions_H
#define OneDFSIDefinitions_H

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// STD
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

// BOOST
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LIFEV
#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/StringUtility.hpp>
#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/algorithm/SolverAmesos.hpp>

namespace ublas = boost::numeric::ublas;

namespace LifeV
{
namespace OneDFSI
{

//#define HAVE_NEUMANN_VISCOELASTIC_BC 1 // Define whether to use homogeneous Neumann/Dirichlet BC for the viscoelastic problem.

/*! @enum Physics Types
 */
enum physicsType_Type
{
    LinearPhysics,        /*!< Use Linear Physics */
    NonLinearPhysics      /*!< Use Non Linear Physics */
};

/*! @enum Flux Types
 */
enum fluxTerm_Type
{
    LinearFlux,        /*!< Use Linear Flux */
    NonLinearFlux      /*!< Use Non Linear Flux */
};

/*! @enum Physics Types
 */
enum sourceTerm_Type
{
    LinearSource,        /*!< Use Linear Source */
    NonLinearSource      /*!< Use Non Linear Source */
};

// Map objects
extern std::map< std::string, physicsType_Type > physicsMap;
extern std::map< std::string, fluxTerm_Type >    fluxMap;
extern std::map< std::string, sourceTerm_Type >  sourceMap;

// Forward class declarations
class OneDFSIModel_Physics;
class OneDFSIModel_Flux;
class OneDFSIModel_Source;
class OneDFSIModel_BCFunction;

enum bcType_Type
{
    W1,         /*!< Riemann variable 1 */
    W2,         /*!< Riemann variable 2 */
    A,          /*!< Area */
    Q,          /*!< Flow rate */
    P,          /*!< Pressure */
    S,          /*!< Normal stress */
    T           /*!< Total normal stress */
};

enum bcSide_Type
{
    left,
    right
};

enum bcLine_Type
{
    first,
    second
};

// ===================================================
// OneDFSIModel Utility Methods
// ===================================================

//! Define the map of the OneDFSIModel objects
inline void
mapsDefinition()
{
    physicsMap["OneD_LinearPhysics"]    = LinearPhysics;
    physicsMap["OneD_NonLinearPhysics"] = NonLinearPhysics;

    fluxMap["OneD_LinearFlux"]          = LinearFlux;
    fluxMap["OneD_NonLinearFlux"]       = NonLinearFlux;

    sourceMap["OneD_LinearSource"]      = LinearSource;
    sourceMap["OneD_NonLinearSource"]   = NonLinearSource;
}

//! Fast pow for the case of exponent 0.5
inline Real
pow05( const Real& base, const Real& exponent )
{
    if ( exponent == 0.5 )
        return std::sqrt( base );
    else
        return std::pow( base, exponent );
}

//! Fast pow for the case of exponent 1.0
inline Real
pow10( const Real& base, const Real& exponent )
{
    if ( exponent == 1.0 )
        return base;
    else
        return std::pow( base, exponent );
}

//! Fast pow for the case of exponent 1.5
inline Real
pow15( const Real& base, const Real& exponent )
{
    if ( exponent == 1.5 )
        return std::sqrt( base ) * base;
    else
        return std::pow( base, exponent );
}

//! Fast pow for the case of exponent 2.0
inline Real
pow20( const Real& base, const Real& exponent )
{
    if ( exponent == 2.0 )
        return base * base;
    else
        return std::pow( base, exponent );
}

//! Fast pow for the case of exponent 3.0
inline Real
pow30( const Real& base, const Real& exponent )
{
    if ( exponent == 3.0 )
        return base * base * base;
    else
        return std::pow( base, exponent );
}

//! Fast pow for the case of exponent 4.0
inline Real
pow40( const Real& base, const Real& exponent )
{
    if ( exponent == 4.0 )
        return base * base * base * base;
    else
        return std::pow( base, exponent );
}

} // OneDFSI namespace
} // LifeV namespace

#endif // OneDFSIDefinitions_H
