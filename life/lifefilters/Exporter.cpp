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

#include <life/lifefilters/Exporter.hpp>

namespace LifeV
{
#if 0
// =================
// Constructor
// =================

template< typename MeshType >
ExporterData<MeshType>::ExporterData( const FieldTypeEnum&   type,
                                      const std::string&     variableName,
                                      const feSpacePtr_Type& feSpacePtr,
                                      const vectorPtr_Type&  vectorPtr,
                                      const UInt&            start,
                                      const FieldRegimeEnum& regime,
                                      const WhereEnum&       where ):
        M_variableName      ( variableName ),
        M_feSpacePtr        ( feSpacePtr ),
        M_storedArrayPtr    ( vectorPtr ),
        M_numDOF            ( feSpacePtr->dim() ),
        M_start             ( start ),
        M_fieldType         ( type ),
        M_regime            ( regime ),
        M_where             ( where )
{}

// ==============
// Operators
// ==============

template< typename MeshType >
Real ExporterData<MeshType>::operator()( const UInt i ) const
{
    return (*M_storedArrayPtr)[i];
}

template< typename MeshType >
Real& ExporterData<MeshType>::operator()( const UInt i )
{
    return (*M_storedArrayPtr)[i];
}

template< typename MeshType >
std::string ExporterData<MeshType>::typeName() const
{
    switch (M_fieldType)
    {
    case ScalarField:
        return "Scalar";
    case VectorField:
        return "Vector";
    }

    return "ERROR string";
}

template< typename MeshType >
UInt ExporterData<MeshType>::typeDim() const
{
    /*switch ( M_fieldType )
    {
    case ScalarField:
        return 1;
    case VectorField:
        return 3;
    }*/

    return M_feSpacePtr->fieldDim();
}

template< typename MeshType >
std::string ExporterData<MeshType>::whereName() const
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
#endif
}
