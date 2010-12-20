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
  @brief exporter and ExporterData classes provide interfaces for post-processing

  @date 11-11-2008
  @author Simone Deparis <simone.deparis.epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>

    Usage: two steps
    <ol>
        <li> first: add the variables using addVariable
        <li> second: call postProcess( time );
    </ol>
*/

#include <life/lifefilters/exporter.hpp>

namespace LifeV
{

// =================
// Constructor
// =================

ExporterData::ExporterData( const ExporterData::Type& type,
                            const std::string& variableName,
                            const vectorPtr_Type& vec,
                            const UInt& start,
                            const UInt& size,
                            const UInt& steady,
                            const ExporterData::Where& where):
        M_variableName  ( variableName ),
        M_vr            ( vec ),
        M_size          ( size ),
        M_start         ( start ),
        M_type          ( type ),
        M_steady        ( steady ),
        M_where         ( where )
{}

// ==============
// Operators
// ==============

Real ExporterData::operator()( const UInt i ) const
{
    return (*M_vr)[i];
}

Real& ExporterData::operator()( const UInt i )
{
    return (*M_vr)[i];
}

std::string ExporterData::typeName() const
{
    switch (M_type)
    {
    case Scalar:
        return "Scalar";
    case Vector:
        return "Vector";
    }

    return "ERROR string";
}

UInt ExporterData::typeDim() const
{
    switch ( M_type )
    {
    case Scalar:
        return 1;
    case Vector:
        return 3;
    }

    return 0;
}

std::string ExporterData::whereName() const
{
    switch (M_where)
    {
    case Node:
        return "Node";
    case Cell:
        return "Cell";
    }

    return "ERROR string";
}

}
