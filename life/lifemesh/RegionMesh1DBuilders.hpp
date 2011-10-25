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
    @brief Contains methods to generate 1D meshes.

    @author -
    @contributor Mauro Perego <mperego@fsu.edu>
    @maintainer -

    @date 10-04-2011
 */

#ifndef REGIONMESH1DBUILDER_HPP
#define REGIONMESH1DBUILDER_HPP 1

#include <life/lifecore/LifeV.hpp>
#include <life/lifemesh/RegionMesh.hpp>
#include <life/lifemesh/MeshChecks.hpp>
#include <fstream>

namespace LifeV
{

//! Build uniform mesh along the x axis.
    /**
     *  @param mesh Reference to the mesh
     *  @param x_l Left end point
     *  @param x_r Right end point
     *  @param numberOfElements Number of elements inside the mesh.
     *
     *  Build 1D uniform mesh along the x axis, extending from x_l to x_r, with numberOfElements elements
     */
template <typename MC>
void uniformMesh1D( RegionMesh<LinearLine, MC>& mesh, const Real& x_l, const Real& x_r, const UInt& numberOfElements )
{
	typedef RegionMesh<LinearLine, MC> mesh_Type;
	ASSERT_PRE( x_r > x_l, "uniformMesh1D problems! The left end point coordinate value should be less than right end point coordinate!");
    ASSERT_PRE( numberOfElements > 0, "uniformMesh1D problems! The number of elements must be positive!");

    mesh.setMaxNumPoints(numberOfElements + 1, true);
    mesh.setNumBPoints(2);

    Real deltax = (x_r - x_l) / numberOfElements;

    typename mesh_Type::point_Type * pp = 0;

    for (UInt it = 0; it < numberOfElements + 1; it++)
    {
        // insert a new Point1D in point list
        pp = &mesh.addPoint( (it == numberOfElements) || ( it == 0) );
        pp->x() = x_l + it*deltax;
        pp->y() = pp->z() = 0.;
        pp->setId(it);
        pp->setLocalId(it);
    }

    mesh.setMaxNumEdges(numberOfElements, true);
    mesh.setNumGlobalVertices( mesh.pointList.size() );
    mesh.setNumVertices(mesh.pointList.size() );

    typename mesh_Type::edge_Type* pe = 0;

    for (UInt it = 0; it < numberOfElements; it++)
    {
        pe = &mesh.addEdge( false );
        pe->setPoint(0, mesh.point(it));
        pe->setPoint(1, mesh.point(it + 1));
        pe->setId(it);
        pe->setLocalId(it);
    }
    mesh.setNumEdges(mesh.edgeList.size() );
    mesh.setMaxNumGlobalEdges(mesh.edgeList.size() );
}



} // Namespace LifeV

#endif /* REGIONMESH1DBUILDER_HPP */
