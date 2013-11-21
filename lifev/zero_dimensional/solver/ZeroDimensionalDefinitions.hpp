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
 *  @brief Zero Dimensional Model Global Definitions
 *  @version alpha (experimental)
 *
 *  @date 06-02-2012
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @mantainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ZeroDimensionalDefinitions_H
#define ZeroDimensionalDefinitions_H


// STD
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>

#ifdef _MSC_VER
#include <iso646.h>
#endif

// BOOST
#include <boost/shared_ptr.hpp>


// LIFEV
#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/util/StringUtility.hpp>
#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/filter/GetPot.hpp>

#define BC_CONSTANT 1000

#define ZERO_DIMENTIONAL_DEFINED_ELEMENTS 6
#define ZERO_DIMENTIONAL_DEFINED_NODES    2

namespace LifeV
{

enum ZeroDimensionalElementType
{
    resistor,
    capacitor,
    inductor,
    diode,
    voltageSource,
    currentSource
};

enum ZeroDimensionalNodeType
{
    knownNode,
    unknownNode
};

enum ZeroDimensionalBCType
{
    Current,
    Voltage
};

} // LifeV namespace

#endif // ZeroDimensionalDefinitions_H
