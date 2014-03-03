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
    @brief MapEpetraData

    @author Antonio Cervone <ant.cervone@gmail.com>
    @maintainer Antonio Cervone <ant.cervone@gmail.com>

    @date 2013-04-19

    This class stores information needed to build a MapEpetra object
 */

#ifndef EPETRAMAPDATA_HPP
#define EPETRAMAPDATA_HPP

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/EnumMapEpetra.hpp>

namespace LifeV
{

//! MapEpetraData
/*!
  The MapEpetraData class stores local list of IDs to build Unique and Repeated maps.
 */
struct MapEpetraData
{
    typedef std::vector<Int> idList_Type;
    idList_Type unique;
    idList_Type repeated;
};

} // namespace LifeV

#endif // EPETRAMAPDATA_HPP
