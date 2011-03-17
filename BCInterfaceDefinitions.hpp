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

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// STL classes
#include <sstream>
#include <string>
#include <vector>

// Boost classes
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LifeV classes
#include <life/lifecore/LifeV.hpp>
#include <life/lifefilters/GetPot.hpp>
#include <life/lifecore/Displayer.hpp>
#include <life/lifecore/Factory.hpp>
#include <life/lifecore/FactorySingleton.hpp>

// 1D BCHandler
#include <lifemc/lifefem/OneDimensionalBCHandler.hpp>

// 3D BCHandler
#include <life/lifefem/BCHandler.hpp>

namespace LifeV
{

enum baseList1D_Type
{
    BCI1DFunction,
    BCI1DFunctionFile,
    BCI1DFunctionSolver,
    BCI1DFunctionFileSolver,
    BCI1DFunctionDefault
};

// Enum objects
enum baseList3D_Type
{
    BCI3DFunction,
    BCI3DFunctionFile,
    BCI3DFunctionSolver,
    BCI3DFunctionFileSolver,
    BCI3DFunctionFSI
};

} // Namespace LifeV

#endif /* BCInterfaceDefinitions_H */
