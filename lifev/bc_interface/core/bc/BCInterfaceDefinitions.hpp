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
 *  @brief File containing the BCInterface definitions
 *
 *  @date 12-11-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceDefinitions_H
#define BCInterfaceDefinitions_H 1


// STL classes
#include <sstream>
#include <string>
#include <vector>

// Boost classes
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>


// LifeV classes
#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>

namespace LifeV
{

enum baseList_Type
{
    BCIFunctionParser,
    BCIFunctionParserFile,
    BCIFunctionParserSolver,
    BCIFunctionParserFileSolver,
    BCIFunctionUserDefined,
    BCIFunctionSolverDefined,
    BCI3DDataInterpolator
};

enum baseContainer_Type
{
    BASEDefault,
    BASEFunction1D,
    BASEFunction3D,
    BASEVector3D,
    BASEVectorInterface3D
};

// Forward class declarations
template< typename BcHandlerType, typename PhysicalSolverType >
class BCInterfaceFactory;

} // Namespace LifeV

#endif /* BCInterfaceDefinitions_H */
