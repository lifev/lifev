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
    @file
    @brief GetPot2Teuchos converter

    @author Umberto Villa <uvilla@emory.edu>
    @contributor
    @maintainer

    @date 04-10-2010

	This utility converts a GetPot data file into a Teuchos XML data file.
	-- boolean variable must be described by "true", "t", "false", "f" (case insensitive)
	-- string are converted to lower case
	-- real variable must contain the character '.' or 'e'
	-- semicolon are ignored
 */

#ifndef FILTER_HPP_
#define FILTER_HPP_

#include<iosfwd>
#include<Teuchos_ParameterList.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifefilters/GetPot.hpp>

/*
 *  *
 */

namespace LifeV
{

//! Type of Variable
/*!
 * Parser Rules:
 * -- IntegerVariable --> only digits and at most one "-"
 * -- RealVariable    --> only digits, at most two of "-" and "+", at exactly one "." or one "e"
 * -- BooleanVariable --> true: true, True, TRUE, t, T (case insensitive)
 *                        false: false, False, FALSE, f, F (case insensitive)
 * -- StringVariable  --> at least one alphabet letter (not interpreted as BooleanVariable)
 * 						  the string will be automatically converted to lower case.
 */
enum GetPotVariable {IntegerVariable, RealVariable, BooleanVariable, StringVariable};

//! Guess the type of the variable contained in the string.
/*!
 * stringValue is modified:
 * 	-- remove semicolon (;) at the end if present
 *  -- convert to lower case
 */
GetPotVariable guessMyType(std::string & stringValue);

//! read the information contained in dataFileName with GetPot and fills a Teuchos Parameter list.
bool fillFromGetPot( const std::string & dataFileName, Teuchos::ParameterList & _pList);


//! convert a string in T object
template <class T>
bool from_string(T& t,
                 const std::string& s)
{
    std::istringstream iss(s);
    return !(iss >> t).fail();
}

}

#endif /* FILTER_HPP_ */
