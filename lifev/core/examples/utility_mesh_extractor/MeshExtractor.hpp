/*
 * MeshExtractor.hpp
 *
 *  Created on: Jan 29, 2012
 *      Author: uvilla
 */

#ifndef MESHEXTRACTOR_HPP_
#define MESHEXTRACTOR_HPP_

#include <lifev/core/mesh/RegionMesh.hpp>

namespace LifeV
{

template<typename geoShape_Type>
RegionMesh< typename RegionMesh<geoShape_Type>::facetShape_Type >*   extractBoundaryMesh (const RegionMesh<geoShape_Type>& mesh3D,
        const UInt& boundaryFaceMarker, const std::list<UInt>& otherBoundaryFaceMarkerList)
{
    //(0) Some useful typedefs
    typedef RegionMesh<geoShape_Type>  mesh_Type;
    typedef typename mesh_Type::facetShape_Type facetShape_Type;
    typedef typename mesh_Type::facet_Type facet_Type;
    typedef typename mesh_Type::ridge_Type ridge_Type;
    typedef RegionMesh<facetShape_Type> boundaryRegionMesh_Type;

    //(1) Extract the list of faces with marker boundaryFaceMarker.
    std::vector<facet_Type const*> boundaryFaces;
    Predicates::EntityMarkerIDInterrogator<facet_Type> predicator (boundaryFaceMarker);
    boundaryFaces = mesh3D.faceList.extractAccordingToPredicate (predicator);

    //(2) Extract the list of points
    std::map<UInt, UInt> vertexIdMap, inverseVertexIdMap; // vertexIdMap[PointIdIn3DMesh] = NewPointIdIn2DMesh
    for ( typename std::vector<facet_Type const*>::iterator it (boundaryFaces.begin() ); it != boundaryFaces.end(); ++it )
    {
        for (UInt i (0); i < boundaryRegionMesh_Type::face_Type::S_numPoints; ++i)
        {
            UInt originalId = (*it)->point (i).id();
            vertexIdMap.insert (std::pair<UInt, UInt> (originalId, vertexIdMap.size() ) );
            inverseVertexIdMap.insert (std::pair<UInt, UInt> (inverseVertexIdMap.size(), originalId) );
        }
    }

    //(3) Extract the list of points shared between the side I want to extract and its neighbor sides
    std::map<UInt, UInt> boundaryPointsOnEdges;
    //boundaryPointsOnEdges[PointIdIn3DMesh] = marker_of_the_other_side_who_shares_the_point
    for (std::list<UInt>::const_iterator it (otherBoundaryFaceMarkerList.begin() ); it != otherBoundaryFaceMarkerList.end(); ++it)
    {
        UInt marker = *it;
        Predicates::EntityMarkerIDInterrogator<facet_Type> predicator (marker);
        std::vector<facet_Type const*> otherboundaryFaces;
        otherboundaryFaces = ( mesh3D.faceList.extractAccordingToPredicate (predicator) );
        std::cout << "Extract " << otherboundaryFaces.size() << " bd faces on Marker" << marker << "\n";
        for ( typename std::vector<facet_Type const*>::iterator it1 (otherboundaryFaces.begin() ); it1 != otherboundaryFaces.end(); ++it1 )
        {
            for (UInt i (0); i < boundaryRegionMesh_Type::face_Type::S_numPoints; ++i)
            {
                UInt originalId = (*it1)->point (i).id();
                if (vertexIdMap.count (originalId) )
                {
                    boundaryPointsOnEdges.insert (std::pair<UInt, UInt> (originalId, marker) );
                }
            }
        }
    }

    //(4) Allocate memory for the new mesh
    boundaryRegionMesh_Type* mesh2D (new boundaryRegionMesh_Type ( mesh3D.comm() ) );
    typename boundaryRegionMesh_Type::point_Type* pp = 0;
    typename boundaryRegionMesh_Type::edge_Type* pe = 0;
    typename boundaryRegionMesh_Type::face_Type* pf = 0;

    UInt __nv (vertexIdMap.size() ), __ne (boundaryPointsOnEdges.size() ), __nt (boundaryFaces.size() );

    debugStream (8000) << "number of vertices: " << __nv << "\n";
    debugStream (8000) << "number of edges: " << __ne << "\n";
    debugStream (8000) << "number of triangles: " << __nt << "\n";

    //(5) Dump the BoundaryPoints coordinates and marker in a std::vector
    std::vector<Real> __x (3 * __nv);
    std::vector<bool> __isonboundary (__nv);
    std::vector<UInt> __whichboundary (__nv);

    // count the number of nodes on the boundary
    UInt __nbv ( 0 );

    // reading vertices
    for ( std::map<UInt, UInt>::iterator it (vertexIdMap.begin() ); it != vertexIdMap.end(); ++it)
    {
        UInt pointId (it->second);
        for (UInt icoor (0); icoor < 3; ++icoor)
        {
            __x[3 * pointId + icoor] = mesh3D.point (it->first).coordinatesArray() [icoor];
        }

        if ( boundaryPointsOnEdges.count (it->first) )
        {
            __whichboundary[ pointId ] = boundaryPointsOnEdges[it->first];
            __isonboundary[ pointId ] = true;
            __nbv                += 1;
        }
        else
        {
            __whichboundary[ pointId ] = 0;
            __isonboundary[ pointId ] = false;
        }
    }

    //(5) Dump the (2D) element connectivity and markers in a std::vector
    std::vector<int> __triangle_nodes ( 3 * __nt );
    std::vector<int> __triangle_label ( __nt );

    MeshElementBareHandler<BareEdge> _be;
    std::pair<BareEdge, bool> _edge;
    UInt counter = 0;
    UInt nTotalEdges (0);
    for ( typename std::vector<facet_Type const*>::iterator it (boundaryFaces.begin() ); it != boundaryFaces.end(); ++it, ++counter )
    {
        for (UInt i (0); i < boundaryRegionMesh_Type::face_Type::S_numPoints; ++i)
        {
            UInt originalId = (*it)->point (i).id();
            __triangle_nodes[ boundaryRegionMesh_Type::face_Type::S_numPoints * counter + i] = vertexIdMap[originalId];
        }

        std::pair<UInt, bool> _check;

        UInt i1 = __triangle_nodes[ 3 * counter ];
        UInt i2 = __triangle_nodes[ 3 * counter + 1 ];
        UInt i3 = __triangle_nodes[ 3 * counter + 2 ];

        _edge                             = makeBareEdge ( i1, i2 );
        _check                            = _be.addIfNotThere ( _edge.first );
        if (!_check.second)
        {
            _be.deleteIfThere ( _edge.first );
        }
        else
        {
            ++nTotalEdges;
        }

        _edge                             = makeBareEdge ( i2, i3 );
        _check                            = _be.addIfNotThere ( _edge.first );
        if (!_check.second)
        {
            _be.deleteIfThere ( _edge.first );
        }
        else
        {
            ++nTotalEdges;
        }

        _edge                             = makeBareEdge ( i3, i1 );
        _check                            = _be.addIfNotThere ( _edge.first );
        if (!_check.second)
        {
            _be.deleteIfThere ( _edge.first );
        }
        else
        {
            ++nTotalEdges;
        }

        __triangle_label[ counter ] = 1;

    }

    debugStream (8000) << "Found " << _be.size() << " boundary edges \n";
    debugStream (8000) << "Found " << nTotalEdges - _be.size() << " internal edges \n";
    __ne = _be.size();

    //(6) Dump the (1D) boundary facet connectivity and markers in a std::vector
    // Here there is some bug, I get too many boundary facets
    std::vector<int> __edge_nodes;
    __edge_nodes.reserve ( 2 * __ne );
    std::vector<int> __edge_label;
    __edge_label.reserve ( __ne );

    for (MeshElementBareHandler<BareEdge>::iterator it (_be.begin() ); it != _be.end(); ++it)
    {
        BareEdge edge (it->first);
        UInt first2dId  = edge.first;
        UInt second2dId = edge.second;
        __edge_nodes.push_back (first2dId);
        __edge_nodes.push_back (second2dId);
        __edge_label.push_back ( boundaryFaceMarker );
    }

    // (7) Set mesh properties
    mesh2D->setMarkerID ( 1 );

    mesh2D->setMaxNumEdges      ( nTotalEdges );
    mesh2D->setMaxNumGlobalEdges ( nTotalEdges );

    mesh2D->setNumEdges         ( nTotalEdges );

    mesh2D->setNumBEdges        ( __ne );

    mesh2D->setMaxNumFaces      ( __nt );
    mesh2D->setMaxNumGlobalFaces ( __nt );

    mesh2D->setNumFaces          ( __nt);

    mesh2D->setMaxNumPoints      ( __nv, true );
    mesh2D->setMaxNumGlobalPoints ( __nv );

    mesh2D->setNumVertices       ( __nv );
    mesh2D->setNumGlobalVertices ( __nv );

    mesh2D->numBVertices()      = __nbv;
    mesh2D->setNumBPoints ( mesh2D->numBVertices() );

    debugStream (8000) << "number of points : " << mesh2D->numPoints() << "\n";
    debugStream (8000) << "number of edges : " << mesh2D->numEdges() << "\n";
    debugStream (8000) << "number of boundary edges : " << mesh2D->numBEdges() << "\n";
    debugStream (8000) << "number of vertices : " << mesh2D->numVertices() << "\n";
    debugStream (8000) << "number of boundary vertices : " << mesh2D->numBVertices() << "\n";

    // (8) Fill in the mesh
    // add points to the mesh
    for ( UInt __i = 0; __i < __nv; ++__i )
    {
        pp = & (mesh2D->addPoint ( __isonboundary[ __i ], false ) );
        pp->setMarkerID ( __whichboundary[ __i ] );
        pp->x() = __x[ 3 * __i ];
        pp->y() = __x[ 3 * __i + 1 ];
        pp->z() = __x[ 3 * __i + 2];
        pp->setId ( __i );
    }

    // add the edges to the mesh
    for ( UInt __i = 0; __i < __ne; ++__i )
    {
        pe = & ( mesh2D->addEdge ( true ) );
        pe->setId ( __i );
        pe->setMarkerID ( markerID_Type ( __edge_label[ __i ] ) );
        pe->setPoint ( 0, mesh2D->point ( __edge_nodes[ 2 * __i ] ) );
        pe->setPoint ( 1, mesh2D->point ( __edge_nodes[ 2 * __i + 1 ] ) );
    }

    // add the triangles to the mesh
    for ( UInt __i = 0; __i < __nt; ++__i )
    {
        pf = & ( mesh2D->addFace (true) );
        pf->setId     ( __i );
        pf->setMarkerID ( markerID_Type ( __triangle_label[ __i ] ) );
        pf->setPoint ( 0, mesh2D->point ( __triangle_nodes[ 3 * __i ] ) );
        pf->setPoint ( 1, mesh2D->point ( __triangle_nodes[ 3 * __i + 1 ] ) );
        pf->setPoint ( 2, mesh2D->point ( __triangle_nodes[ 3 * __i + 2 ] ) );
    }


    //(9) Update all the facets
    mesh2D->updateElementFacets (true);


    return mesh2D;
}
}


#endif /* MESHEXTRACTOR_HPP_ */
