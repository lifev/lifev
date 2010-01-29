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
 *  @brief Parser Definitions
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 29-01-2010
 */

#ifndef Parser_Definitions_H
#define Parser_Definitions_H 1

// STD Classes
#include <map>
#include <iomanip>
#include <string>

// LifeV config
#include <lifeconfig.h>
#include <life/lifecore/life.hpp>

// BOOST Classes
#include <boost/algorithm/string.hpp>
#ifdef HAVE_BOOST_SPIRIT_CLASSIC
    #ifdef HAVE_BOOST_SPIRIT_CLASSIC_TEMP
        #include <boost/spirit/include/classic.hpp>
        #include <boost/spirit/include/phoenix1_binders.hpp>

        namespace spirit = boost::spirit::classic;
    #else
        #include <boost/spirit.hpp>
        #include <boost/spirit/phoenix/binders.hpp>

        namespace spirit = boost::spirit;
    #endif
#endif

namespace LifeV {

typedef std::map< std::string, Real >            Variables_Type;
typedef std::vector< Real >                      Results_Type;
typedef std::vector< std::string >               Strings_Type;

} // Namespace LifeV

#endif /* Parser_Definitions_H */
