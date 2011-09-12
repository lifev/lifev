/** @file Contains utility for importing meshes
 * BareMesh.hpp
 *
 *      @date 16 Aug 2011
 *      @author: luca formaggia
 */
#include <life/lifecore/LifeV.hpp>
#include <map>
#include <vector>
#include <life/lifemesh/ElementShapes.hpp>
#include <life/lifearray/ArraySimple.hpp>

#ifndef BAREMESH_HPP_
#define BAREMESH_HPP_

namespace LifeV
{

//! A struct for a bare mesh
/**
 * A very simple struct which stores an mesh as read from a file, ready to be imported in
 * a regionmesh
 * All SimpleArray have the first dimension the "shortest" one
 */
template <typename GeoShape>
struct RegionMeshBare{
    UInt nDimensions;
    ReferenceShapes refShape;
    ReferenceShapes bRefShape;
    ArraySimple<Real> points;
    std::vector<ID> pointsMarkers;
    ArraySimple<UInt> edges;
    std::vector<ID> edgesMarkers;
    UInt numBoundaryFaces;
    ArraySimple<UInt> faces;
    std::vector<ID> facesMarkers;
    ArraySimple<UInt> elements;
    std::vector<ID> elementsMarkers;

    void clear();
};

template <typename GeoShape>
void RegionMeshBare<GeoShape>::clear()
{
    clearVector ( points );
    clearVector ( pointsMarkers );
    clearVector ( edges );
    clearVector ( edgesMarkers );
    clearVector ( faces );
    clearVector ( facesMarkers );
    clearVector ( elements );
    clearVector ( elementsMarkers );
}

}

#endif /* BAREMESH_HPP_ */
