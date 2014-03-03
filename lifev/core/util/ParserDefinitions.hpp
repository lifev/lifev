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
 *  @brief File containing the Parser definitions
 *
 *  @date 29-01-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef Parser_Definitions_H
#define Parser_Definitions_H 1


// STD Classes
#include <iomanip>
#include <algorithm>

// BOOST Classes
#include <boost/algorithm/string.hpp>

#ifdef HAVE_BOOST_SPIRIT_QI
#ifdef ENABLE_SPIRIT_PARSER

// BOOST Spirit Classes
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

// BOOST Spirit Namespaces
namespace spirit  = boost::spirit;
namespace qi      = spirit::qi;
namespace ascii   = spirit::ascii;
namespace phoenix = boost::phoenix;

#endif /* ENABLE_SPIRIT_PARSER */
#endif /* HAVE_BOOST_SPIRIT_QI */


// LifeV classes
#include <lifev/core/LifeV.hpp>

#endif /* Parser_Definitions_H */
