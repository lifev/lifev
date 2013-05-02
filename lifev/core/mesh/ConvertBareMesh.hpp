/********************************************************************************
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
********************************************************************************/

/**
 * @file   ConvertBareMesh.hpp
 * @brief  Convert a BareMesh into a RegionMesh.
 * @author Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 * @date   10/2012
**/

#ifndef CONVERT_BARE_MESH_HPP__
#define CONVERT_BARE_MESH_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/MeshElementBare.hpp>
#include <lifev/core/mesh/MeshChecks.hpp>
#include <lifev/core/mesh/BareMesh.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

namespace LifeV
{

namespace
{

template <int dim> struct convert_spec;

template <class S, class MC>
struct converter
{

    typedef BareMesh<S>        bmesh_t;
    typedef RegionMesh<S, MC>  mesh_t;

    typedef typename mesh_t::point_Type    point_t;
    typedef typename mesh_t::ridge_Type    ridge_t;
    typedef typename mesh_t::facet_Type    facet_t;
    typedef typename mesh_t::element_Type  element_t;

    static bool convert (bmesh_t& baremesh, mesh_t& mesh, bool verbose = false)
    {
        // make sure global id are stored when the mesh is unpartitioned
        if ( !baremesh.isPartitioned )
        {
            if ( baremesh.pointIDs.size() == 0 )
            {
                baremesh.pointIDs.resize ( baremesh.points.numberOfColumns() );
                for ( UInt i = 0; i < baremesh.points.numberOfColumns(); ++i )
                {
                    baremesh.pointIDs[i] = i;
                }
            }
            if ( baremesh.ridgeIDs.size() == 0 )
            {
                baremesh.ridgeIDs.resize ( baremesh.ridges.numberOfColumns() );
                for ( UInt i = 0; i < baremesh.ridges.numberOfColumns(); ++i )
                {
                    baremesh.ridgeIDs[i] = i;
                }
            }
            if ( baremesh.facetIDs.size() == 0 )
            {
                baremesh.facetIDs.resize ( baremesh.facets.numberOfColumns() );
                for ( UInt i = 0; i < baremesh.facets.numberOfColumns(); ++i )
                {
                    baremesh.facetIDs[i] = i;
                }
            }
            if ( baremesh.elementIDs.size() == 0 )
            {
                baremesh.elementIDs.resize ( baremesh.elements.numberOfColumns() );
                for ( UInt i = 0; i < baremesh.elements.numberOfColumns(); ++i )
                {
                    baremesh.elementIDs[i] = i;
                }
            }
        }

        // ============
        // begin.POINTS
        // ============
        {
            UInt numberPoints = baremesh.points.numberOfColumns();

            mesh.setIsPartitioned (baremesh.isPartitioned);
            // Update mesh containers
            mesh.setMaxNumPoints (numberPoints, true);
            mesh.setMaxNumGlobalPoints (numberPoints);
            mesh.setNumVertices ( baremesh.numVertices );
            mesh.setNumGlobalVertices ( baremesh.numVertices );
            mesh.setNumBPoints ( baremesh.numBoundaryPoints );
            mesh.setNumBVertices ( baremesh.numBoundaryVertices );
            // Add the points
            for (UInt i = 0; i < numberPoints; ++i)
            {
                // Add a generic point (vertex selection is made later)
                point_t& p = mesh.addPoint (false, false);
                // Coordinates
                p.x() = baremesh.points (0, i);
                p.y() = baremesh.points (1, i);
                p.z() = baremesh.points (2, i);
                // global Id
                p.setId (baremesh.pointIDs[i]);
                // Marker
                p.setMarkerID (baremesh.pointMarkers[i]);
            }
        }
        if (verbose)
            std::cout << "[INFO:convertBareMesh] Added "
                      << mesh.numPoints() << " marked points." << std::endl;

        // ==============
        // begin.ELEMENTS
        // ==============
        {
            // Update containers
            UInt numberElements = baremesh.elements.numberOfColumns();
            mesh.setMaxNumGlobalElements (numberElements);
            mesh.setMaxNumElements (numberElements, true);
            mesh.setNumElements (numberElements);
            for (UInt i = 0; i < numberElements; ++i)
            {
                // Add the element
                element_t& e = mesh.addElement();
                e.setId (baremesh.elementIDs[i]);
                e.setMarkerID (baremesh.elementMarkers[i]);
                // Points
                for (UInt j = 0; j < element_t::S_numPoints; ++j)
                {
                    UInt pid = baremesh.elements (j, i);
                    e.setPoint (j, mesh.point (pid) );
                    // Mark vertices
                    if (j < element_t::S_numVertices)
                    {
                        mesh.pointList (pid).setFlag (LifeV::EntityFlags::VERTEX);
                    }
                }
            }
            // Update number of vertices
            UInt numVertices = mesh.pointList.countElementsWithFlag (
                                   LifeV::EntityFlags::VERTEX,
                                   &LifeV::Flag::testAllSet);
            mesh.setNumVertices (numVertices);
        }
        if (verbose)
            std::cout << "[INFO:convertBareMesh] Added "
                      << mesh.numElements() << " marked elements." << "\n"
                      << "[INFO:convertBareMesh] Found "
                      << mesh.numVertices() << " vertices." << std::endl;

        // Set the region marker
        mesh.setMarkerID (baremesh.regionMarkerID);
        mesh.setId (baremesh.regionMarkerID);

        // dimension specific entities
        convert_spec<S::S_nDimensions>::add_facets (baremesh, mesh, verbose);
        convert_spec<S::S_nDimensions>::add_ridges (baremesh, mesh, verbose);

        // discard baremesh
        baremesh.clear();

        return convert_spec<S::S_nDimensions>::check (mesh, verbose);

    }

};

// Generic conversion (say 3d)
template <int dim>
struct convert_spec
{

    template <class S, class MC>
    static void add_facets (BareMesh<S>& baremesh, RegionMesh<S, MC>& mesh, bool verbose)
    {
        // ============
        // begin.FACETS
        // ============
        typedef RegionMesh<S, MC>  mesh_t;
        typedef typename mesh_t::facet_Type facet_t;
        if (verbose)
            std::cout << "[INFO:convertBareMesh] Updating facets"
                      << std::endl;
        // Update containers
        UInt numberFacets = baremesh.facets.numberOfColumns();
        mesh.setMaxNumFacets (numberFacets);
        mesh.setMaxNumGlobalFacets (numberFacets);
        mesh.setNumFacets (numberFacets);
        mesh.setNumBoundaryFacets (0);
        for (UInt i = 0; i < numberFacets; ++i)
        {
            // Add the facet
            facet_t& f = mesh.addFacet (false);
            f.setId (baremesh.facetIDs[i]);
            f.setMarkerID (baremesh.facetMarkers[i]);
            // Points
            for (UInt j = 0; j < facet_t::S_numPoints; ++j)
            {
                UInt pid = baremesh.facets (j, i);
                f.setPoint (j, mesh.point (pid) );
            }
        }
    }

    template <class S, class MC>
    static void add_ridges (BareMesh<S>& baremesh, RegionMesh<S, MC>& mesh, bool verbose)
    {
        // ============
        // begin.RIDGES
        // ============
        typedef RegionMesh<S, MC>  mesh_t;
        typedef typename mesh_t::ridge_Type ridge_t;
        if (verbose)
            std::cout << "[INFO:convertBareMesh] Updating ridges"
                      << std::endl;
        // Update containers
        UInt numberRidges = baremesh.ridges.numberOfColumns();
        mesh.setMaxNumRidges (numberRidges);
        mesh.setMaxNumGlobalRidges (numberRidges);
        mesh.setNumRidges (numberRidges);
        mesh.setNumBoundaryRidges (0);
        for (UInt i = 0; i < numberRidges; ++i)
        {
            // Add the facet
            ridge_t& r = mesh.addRidge (false);
            r.setId (baremesh.ridgeIDs[i]);
            r.setMarkerID (baremesh.ridgeMarkers[i]);
            // Points
            for (UInt j = 0; j < ridge_t::S_numPoints; ++j)
            {
                UInt pid = baremesh.ridges (j, i);
                r.setPoint (j, mesh.point (pid) );
            }
        }
    }

    template <class S, class MC>
    static bool check (RegionMesh<S, MC>& mesh, bool verbose)
    {
        Switch sw;
        std::stringstream discardedLog;
        std::ostream& oStr = verbose ? std::cout : discardedLog;
        return checkMesh3D (mesh, sw, true, verbose, oStr, std::cerr, oStr);
    }

};

// 2d mesh version
template <>
struct convert_spec<2>
{
    template <class S, class MC>
    static void add_facets (BareMesh<S>& baremesh, RegionMesh<S, MC>& mesh, bool verbose)
    {
        convert_spec<3>::template add_facets (baremesh, mesh, verbose);
    }

    template <class S, class MC>
    static void add_ridges (BareMesh<S>&, RegionMesh<S, MC>&, bool)
    {}

    template <class S, class MC>
    static bool check (RegionMesh<S, MC>&, bool)
    {
        return true;
    }

};

// 1d mesh version
template <>
struct convert_spec<1>
{
    template <class S, class MC>
    static void add_facets (BareMesh<S>&, RegionMesh<S, MC>&, bool)
    {}

    template <class S, class MC>
    static void add_ridges (BareMesh<S>&, RegionMesh<S, MC>&, bool)
    {}

    template <class S, class MC>
    static bool check (RegionMesh<S, MC>&, bool)
    {
        return true;
    }

};

}

//! convertBareMesh - convert a previously read BareMesh in a RegionMesh object
/*!
  Starting from a BareMesh, this routine generates a fully compliant RegionMesh object

  @param bareMesh, the bare mesh data structure in input.
  @param mesh, the mesh data structure to fill in.
  @param verbose, setting it as true, the output is verbose (the default is false).
  @return true if everything went fine, false otherwise.
*/
template <typename GeoShapeType, typename MCType>
bool
convertBareMesh ( BareMesh<GeoShapeType>& bareMesh,
                  RegionMesh<GeoShapeType, MCType>&  mesh,
                  bool                       verbose = false)
{
    return converter<GeoShapeType, MCType>::convert (bareMesh, mesh, verbose);
}// Function convertBareMesh

} // Namespace LifeV

#endif // CONVERT_BARE_MESH_HPP__
