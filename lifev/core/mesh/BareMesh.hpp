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
    @file BareMesh.hpp
    @brief Contains utility for importing meshes

    @date 16 Aug 2011
    @author: luca formaggia
*/

#ifndef BAREMESH_HPP_
#define BAREMESH_HPP_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/ElementShapes.hpp>
#include <lifev/core/array/ArraySimple.hpp>

namespace LifeV
{

//! A struct for a bare mesh
/**
 * A very simple struct which stores an mesh as read from a file, ready to be imported in
 * a regionmesh
 * All SimpleArray have the first dimension the "shortest" one
 */
template <typename GeoShapeType>
struct BareMesh
{
    UInt nDimensions;
    ID regionMarkerID;
    ReferenceShapes refShape;
    ReferenceShapes bRefShape;
    UInt numBoundaryPoints;
    UInt numVertices;
    UInt numBoundaryVertices;
    bool isPartitioned;
    ArraySimple<Real> points;
    std::vector<ID> pointMarkers;
    std::vector<ID> pointIDs;
    ArraySimple<UInt> ridges;
    std::vector<ID> ridgeMarkers;
    std::vector<ID> ridgeIDs;
    UInt numBoundaryFacets;
    ArraySimple<UInt> facets;
    std::vector<ID> facetMarkers;
    std::vector<ID> facetIDs;
    ArraySimple<UInt> elements;
    std::vector<ID> elementMarkers;
    std::vector<ID> elementIDs;

    BareMesh();
    void clear();
};

template <typename GeoShapeType>
BareMesh<GeoShapeType>::BareMesh() : isPartitioned (false) {}

template <typename GeoShapeType>
void BareMesh<GeoShapeType>::clear()
{
    clearVector ( points );
    clearVector ( pointMarkers );
    clearVector ( pointIDs );
    clearVector ( ridges );
    clearVector ( ridgeMarkers );
    clearVector ( ridgeIDs );
    clearVector ( facets );
    clearVector ( facetMarkers );
    clearVector ( facetIDs );
    clearVector ( elements );
    clearVector ( elementMarkers );
    clearVector ( elementIDs );
}

}

#endif /* BAREMESH_HPP_ */
