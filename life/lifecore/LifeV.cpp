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
    @brief definition of nDimensions

    @author Simone Deparis <simone.deparis@epfl.ch>

    @date 2009-12-09

 */

#include <life/lifecore/LifeV.hpp>

namespace LifeV
{
const UInt nDimensions(NDIM);

// flag testers
static bool flagTestAllSet( flag_Type const & inputFlag, flag_Type const & refFlag )
{
  return ( inputFlag  & refFlag ) == refFlag;
  // returns true if all byte-flags common set in refFlag are also set in inputFlag
}

static bool flagTestOneSet( flag_Type const & inputFlag, flag_Type const & refFlag )
{
 return inputFlag  & refFlag;
  // returns true if at least one flag set in refFlag is set in inputFlag
}

} //end namespace LifeV

