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
 *  @brief ExporterData
 *
 *  @author Simone Deparis <simone.deparis@epfl.ch>
 *  @date 11-11-2008
 */

#include <life/lifefilters/exporter.hpp>

namespace LifeV {

ExporterData::ExporterData( const  ExporterData::Type type,
                            const std::string variableName,
                            const vector_ptrtype& vr,
                            UInt start,
                            UInt size,
                            UInt steady ):
    M_variableName  ( variableName ),
    M_vr            ( vr ),
    M_size          ( size ),
    M_start         ( start ),
    M_type          ( type ),
    M_steady        ( steady )
{}

const std::string&
ExporterData::variableName() const
{
    return M_variableName;
}

Real
ExporterData::operator()( const UInt i ) const
{
    return (*M_vr)[i];
}

Real&
ExporterData::operator()( const UInt i )
{
    return (*M_vr)[i];
}

const UInt&
ExporterData::size() const
{
    return M_size;
}

const ExporterData::Type&
ExporterData::type() const
{
    return M_type;
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

UInt
ExporterData::typeDim() const
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

}
