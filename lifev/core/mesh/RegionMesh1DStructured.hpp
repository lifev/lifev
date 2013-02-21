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

    @author Antonio Cervone <ant.cervone@gmail.com>
    @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
    @contributor Mauro Perego <mperego@fsu.edu>
    @maintainer -

    @date 15-10-2012
 */

#ifndef REGIONMESH1DSTRUCTURED_HPP
#define REGIONMESH1DSTRUCTURED_HPP 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

namespace LifeV
{

// Labels for the structured 1D mesh
namespace Structured1DLabel
{
//! Label for the internal entities
const markerID_Type INTERNAL = 0;

//! Label for the left boundary corner
const markerID_Type LEFT = 1;

//! Label for the right boundary corner
const markerID_Type RIGHT = 2;
} // namespace Structured1DLabel

//! Build uniform mesh along the x axis.
/*!
    @tparam mesh Reference to the mesh.
    @param regionFlag flag for the mesh.
    @param numberOfElements Number of elements inside the mesh.
    @param verbose Output verbosity.
    @param lenght Length of the mesh.
    @param origin Origin of the mesh.

    Build 1D uniform mesh along the x axis, extending from origin to origin + length,
    with numberOfElements elements.
*/
template <typename MeshType>
void regularMesh1D ( MeshType& mesh,
                     markerID_Type regionFlag,
                     const UInt& numberOfElements,
                     bool verbose = false,
                     const Real& length = 1.,
                     const Real& origin = 0. )
{
    typedef MeshType mesh_Type;

    if ( verbose && mesh.comm()->MyPID() == 0 )
    {
        std::cout << "Building 1d mesh" << std::endl;
    }

    mesh.setMaxNumPoints (numberOfElements + 1, true);
    mesh.setNumBPoints (2);
    mesh.setMarkerID (regionFlag);

    Real deltax = (length - origin) / numberOfElements;

    typename mesh_Type::point_Type* pp = 0;

    for (UInt it = 0; it < numberOfElements + 1; it++)
    {
        bool isBoundary = (it == numberOfElements) || ( it == 0);

        // insert a new Point1D in point list
        pp = &mesh.addPoint ( isBoundary, false );
        pp->x() = origin + it * deltax;
        pp->y() = 0.;
        pp->z() = 0.;
        pp->setId (it);

        if ( it == 0 )
        {
            pp->firstAdjacentElementIdentity() = 0;
            pp->firstAdjacentElementPosition() = 0;

            pp->secondAdjacentElementIdentity() = NotAnId;
            pp->secondAdjacentElementPosition() = NotAnId;

            pp->setMarkerID ( Structured1DLabel::LEFT );
        }
        else if ( it == numberOfElements )
        {
            pp->firstAdjacentElementIdentity() = it - 1;
            pp->firstAdjacentElementPosition() = 1;

            pp->secondAdjacentElementIdentity() = NotAnId;
            pp->secondAdjacentElementPosition() = NotAnId;

            pp->setMarkerID ( Structured1DLabel::RIGHT );
        }
        else
        {
            pp->firstAdjacentElementIdentity() = it;
            pp->firstAdjacentElementPosition() = 0;

            pp->secondAdjacentElementIdentity() = it - 1;
            pp->secondAdjacentElementPosition() = 1;

            pp->setMarkerID ( Structured1DLabel::INTERNAL );
        }


    }

    mesh.setNumGlobalVertices ( mesh.pointList.size() );
    mesh.setNumVertices (mesh.pointList.size() );
    mesh.setMaxNumPoints      ( mesh.pointList.size(), true );
    mesh.setMaxNumGlobalPoints ( mesh.pointList.size() );
    mesh.numBVertices() = 2;
    mesh.setNumBPoints ( mesh.numBVertices() );

    mesh.setLinkSwitch ( "FACETS_HAVE_ADIACENCY" );
    mesh.setLinkSwitch ( "HAS_ALL_FACETS" );

    mesh.setMaxNumEdges (numberOfElements, true);

    typename mesh_Type::edge_Type* pe = 0;

    for (UInt it = 0; it < numberOfElements; it++)
    {
        pe = &mesh.addEdge ( false );
        pe->setPoint (0, mesh.point (it) );
        pe->setPoint (1, mesh.point (it + 1) );
        pe->setId (it);
        pe->setMarkerID ( Structured1DLabel::INTERNAL );
    }
    mesh.setNumEdges (mesh.edgeList.size() );
    mesh.setMaxNumGlobalEdges (mesh.edgeList.size() );

    mesh.updateElementFacets ( true, false, mesh.pointList.size() );
} // regularMesh1D

} // Namespace LifeV

#endif /* REGIONMESH1DSTRUCTURED_HPP */
