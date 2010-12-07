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

#ifndef BCInterface_Definitions_H
#define BCInterface_Definitions_H 1

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
#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/displayer.hpp>
#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>

#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/bcCond.hpp>
#include <life/lifefem/bcFunction.hpp>
#include <life/lifefem/bcVector.hpp>

namespace LifeV
{

// Enum objects
enum BCInterface_BaseList
{
    BCInterface_function,
    BCInterface_functionFile,
    BCInterface_OPERfunction,
    BCInterface_OPERfunctionFile,
    BCInterface_OPERFSI
};

// Type definitions
typedef EntityFlag        BCFlag;
typedef std::vector< ID > BCComV;

} // Namespace LifeV

#endif /* BCInterface_Definitions_H */
