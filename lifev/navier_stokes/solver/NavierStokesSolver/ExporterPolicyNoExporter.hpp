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
    @file ExporterPolicyNoExporter class
    @brief This class contains is used when no exporter is required

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 04-12-2012
 */

#ifndef EXPORTERPOLICYNOEXPORTER_HPP
#define EXPORTERPOLICYNOEXPORTER_HPP

#include <boost/shared_ptr.hpp>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>


namespace LifeV
{

struct ExporterPolicyNoExporter
{

    typedef VectorEpetra                             vector_Type;
    typedef boost::shared_ptr<VectorEpetra>          vectorPtr_Type;

    static void initExporter ( Teuchos::ParameterList& /*list*/,
                               vectorPtr_Type /*solution*/ ) {}

    static void exportSolution () {}

    static void finalizeExporter() {}

};

} // namespace LifeV

#endif /* EXPORTERPOLICYNOEXPORTER_HPP */
