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
 *  @brief BCInterface Definitions
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 12-11-2009
 */

#ifndef BCInterface_Definitions_H
#define BCInterface_Definitions_H 1

//#define DEBUG 1;

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

// Boost classes
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

// STL classes
#include <sstream>
#include <string>
#include <vector>

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
typedef std::string       BCName;
typedef EntityFlag        BCFlag;
typedef std::vector< ID > BCComV;

} // Namespace LifeV

#endif /* BCInterface_Definitions_H */
