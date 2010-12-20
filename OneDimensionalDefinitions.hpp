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
 *  @mantainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ONEDIMENSIONALMODEL_DEFINITIONS_H
#define ONEDIMENSIONALMODEL_DEFINITIONS_H

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
#include <life/lifecore/life.hpp>
#include <life/lifecore/util_string.hpp>
#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>
#include <life/lifearray/tab.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifealg/SolverAmesos.hpp>

namespace ublas = boost::numeric::ublas;

namespace LifeV
{
namespace OneDimensional
{

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
class OneDimensionalModel_Physics;
class OneDimensionalModel_Flux;
class OneDimensionalModel_Source;
class OneDimensionalModel_BCFunction;

enum bcType_Type
{
    W1,
    W2,
    A,
    Q,
    P
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
// OneDimensionalModel Utility Methods
// ===================================================

//! Define the map of the OneDimensionalModel objects
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

} // OneDimensional namespace
} // LifeV namespace

#endif // ONEDIMENSIONALMODEL_DEFINITIONS_H
