//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

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
 *  @brief One Dimensional Model Global Definitions
 *
 *  @version 1.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 15-04-2010
 */

#ifndef ONEDIMENSIONALMODEL_DEFINITIONS_H
#define ONEDIMENSIONALMODEL_DEFINITIONS_H

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

/*! @enum Physics Types
 */
enum OneDimensionalModel_PhysicsTypes
{
    OneD_LinearPhysics,        /*!< Use Linear Physics */
    OneD_NonLinearPhysics      /*!< Use Non Linear Physics */
};

/*! @enum Flux Types
 */
enum OneDimensionalModel_FluxTypes
{
    OneD_LinearFlux,        /*!< Use Linear Flux */
    OneD_NonLinearFlux      /*!< Use Non Linear Flux */
};

/*! @enum Physics Types
 */
enum OneDimensionalModel_SourceTypes
{
    OneD_LinearSource,        /*!< Use Linear Source */
    OneD_NonLinearSource      /*!< Use Non Linear Source */
};

// Map objects
extern std::map< std::string, OneDimensionalModel_PhysicsTypes > OneDimensionalModel_PhysicsMap;
extern std::map< std::string, OneDimensionalModel_FluxTypes >    OneDimensionalModel_FluxMap;
extern std::map< std::string, OneDimensionalModel_SourceTypes >  OneDimensionalModel_SourceMap;

// Forward class declarations
class OneDimensionalModel_Physics;
class OneDimensionalModel_Flux;
class OneDimensionalModel_Source;
class OneDimensionalModel_BCFunction;

// Type definitions
typedef singleton< factory< OneDimensionalModel_Physics,
OneDimensionalModel_PhysicsTypes > > Factory_OneDimensionalModel_Physics;
typedef singleton< factory< OneDimensionalModel_Flux,
OneDimensionalModel_FluxTypes > >    Factory_OneDimensionalModel_Flux;
typedef singleton< factory< OneDimensionalModel_Source,
OneDimensionalModel_SourceTypes > >  Factory_OneDimensionalModel_Source;

//typedef singleton< factoryClone< OneDimensionalModel_BCFunction > > FactoryClone_OneDimensionalModel_BCFunction;

typedef boost::array< Real, 2 >                 Container2D_Type;

// ScalVec SHOULD BE REPLACED EVERYWHERE BY EPETRAVECTOR FOR PARALLEL COMPUTATION
typedef ublas::vector< Real >                   ScalVec;

enum OneD_BC
{
    OneD_W1,
    OneD_W2,
    OneD_A,
    OneD_Q,
    OneD_P
};

enum OneD_BCSide
{
    OneD_left,
    OneD_right
};

enum OneD_BCLine
{
    OneD_first,
    OneD_second
};

// ===================================================
// OneDimensionalModel Utility Methods
// ===================================================

//! Define the map of the OneDimensionalModel objects
inline void
OneDimensionalModel_MapsDefinition()
{
    OneDimensionalModel_PhysicsMap["OneD_LinearPhysics"]    = OneD_LinearPhysics;
    OneDimensionalModel_PhysicsMap["OneD_NonLinearPhysics"] = OneD_NonLinearPhysics;

    OneDimensionalModel_FluxMap["OneD_LinearFlux"]          = OneD_LinearFlux;
    OneDimensionalModel_FluxMap["OneD_NonLinearFlux"]       = OneD_NonLinearFlux;

    OneDimensionalModel_SourceMap["OneD_LinearSource"]      = OneD_LinearSource;
    OneDimensionalModel_SourceMap["OneD_NonLinearSource"]   = OneD_NonLinearSource;
}

//! Scalar product between 2D vectors
inline Real
dot( const Container2D_Type& vector1,
     const Container2D_Type& vector2 )
{
    ASSERT_PRE( vector1.size() == 2 && vector2.size() == 2, "dot works only for 2D vectors" );

    return vector1[0]*vector2[0] + vector1[1]*vector2[1];
}

}

#endif // ONEDIMENSIONALMODEL_DEFINITIONS_H
