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
    @brief Contains the basic element shapes, to be used by Geometric
           and Finite Elements

    @author Luca Formaggia
    @contributor Zhen Wang <zwang26@emory.edu>
    @contributor Tiziano Passerini <tiziano@mathcs.emory.edu>

    These classes provide the most basic definitions for the geometric
    elements. The variables names are self-explanatory.

    The Basis Geometric Shapes (<tt>GeoShape</tt>) are derived from a set of
    classes called Basis Reference Shapes (<tt>BasRefSha</tt>), which contains
    the very basic information about the reference element geometries
    supported by the library (reference element geometry=geometry of the
    finite element in the reference space.)

    A Basis Geometric Shape contains (when relevant) the following methods
    <ol>
    <li> <tt>static ID edgeToPoint(ID const jEdge, ID const iPoint);</tt>
         which returns the local ID of the i-th point of the j-th edge of the GeoShape

    <li> <tt>static ID faceToPoint(ID const jFace, ID const iPoint);</tt>
         which returns the local ID of the i-th point of the j-th face of the GeoShape

    <li> <tt>static pair<ID,bool> faceToEdge(ID const jFace, ID const iEdge);</tt>
         which returns the local numbering of the i-th edge on the j-th face.
         It returns also if the orientation of the edge on the face is consistent
         with that of the same edge on the element
    </ol>

    @note The methods edge to point-ID (EtoP) and Face to point-ID (FtoP)
          return the local id number of points on faces and edges (when relevant)
    @note We follow the convention of indicating THE VERTICES FIRST in the list
          of dofs
 */

#ifndef ELEMENTSHAPES_H
#define ELEMENTSHAPES_H 1

#include <utility>
#include <lifev/core/LifeV.hpp>

namespace LifeV
{
//! A utility to invert point numbering on a GeoShape
/*!
 *  It must be specialised for the specific GeoShape since the inverse ordering
 *  depends on how points are numbered in the actual GeoShape.
 *  This utility is meant to be used only by procedures that build a mesh, since it operates on
 *  basic mesh structures. It can be dangerous to use, for instance, after a full mesh has been set up.
 *  It is useful to invert faces or edges which are incorrectly oriented or to fix a mesh produced by a mesher
 *  which uses a different orientation convention.
 *
 *  @param pointId Elemental local id of a point of the GeoShape
 *  @return the (local) ID of the corresponding point in the reversed GeoShape
 */
template <typename GeoShapeType>
inline ID reversePoint ( ID const& pointId );

//                   *********** BASIS REFERENCE SHAPES ****************

/*! @enum ReferenceShapes
    Lists of the geometries of the finite element in the reference space,
    supported by the library

enum ReferenceShapes
{
    NONE, POINT, LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA
};
 */
// deprecated name
enum ReferenceShapes
{
    NONE, POINT, LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA
};


//!
/*!
    @param shape a shape identifier
    @return the geometric dimension of the shape
    @sa ReferenceShapes
 */
UInt shapeDimension (const ReferenceShapes& shape);


/*! @enum ReferenceGeometry
    Lists of the geometric items used to build the shapes.

enum ReferenceGeometry
{
    VERTEX = 0, EDGE = 1, FACE = 2, VOLUME = 3
};
 */
// deprecated name
enum ReferenceGeometry
{
    VERTEX = 0, EDGE = 1, FACE = 2, VOLUME = 3
};


//! @defgroup BasRefSha   Basis Reference Shapes

//! @ingroup BasRefSha
class nullShape
{
public:
    typedef nullShape GeoBShape;
};


//! dummy class for selecting correct function specializations based on geometry dimensions (1,2,3).
template<int geoDimensions>
class GeoDim {};


//! @ingroup BasRefSha
class Point
{
public:
    static const ReferenceShapes S_shape      = POINT; //!< Identify the shape
    static const ReferenceGeometry S_geometry = VERTEX;//!< Identify the geometric entity
    static const UInt S_nDimensions           = 0;     //!< Dimensionality
    static const UInt S_numFaces              = 0;     //!< Number of faces
    static const UInt S_numEdges              = 0;     //!< Number of edges
    static const UInt S_numVertices           = 1;     //!< Number of vertices.
    static const UInt S_numFacets             = 0;     //!< Number of facets
    static const UInt S_numRidges             = 0;     //!< Number of ridges
    static const UInt S_numPeaks              = 0;     //!< Number of peaks

    //! @return the local ID of the j-th point of the i-th edge
    static ID edgeToPoint ( ID const& /*iEdge*/, ID const& /*jPoint*/ )
    {
        ERROR_MSG ( "edgeToPoint not implemented for geo Point elements." );
        return static_cast<ID> (0);
    }
};


//! @ingroup BasRefSha
class Line
{
public:
    static const ReferenceShapes S_shape = LINE;     //!< Identify the shape
    static const ReferenceGeometry S_geometry = EDGE;//!< Identify the geometric entity
    static const UInt S_nDimensions = 1;             //!< Dimensionality
    static const UInt S_numFaces = 0;                //!< Number of faces
    static const UInt S_numEdges = 1;                //!< Number of edges
    static const UInt S_numVertices = 2;             //!< Number of vertices.
    static const UInt S_numFacets = S_numVertices;   //!< Number of facets
    static const UInt S_numRidges = 0;               //!< Number of ridges
    static const UInt S_numPeaks = 0;                //!< Number of peaks

    //! @return the local ID of the j-th point of the i-th edge
    static ID edgeToPoint ( ID const& /*iEdge*/, ID const& jPoint )
    {
        return jPoint;
    }

    //! @return the local ID of the j-th point of the i-th facet
    static ID facetToPoint ( ID const& iFacet, ID const& /*jPoint*/ )
    {
        return iFacet;
    }

    static ID facetToRidge ( ID const& /*iFacet*/, ID const& /*jRidge*/ )
    {
        return NotAnId;
    }

    static ID facetToPeak ( ID const& /*iFacet*/, ID const& /*jPeak*/ )
    {
        return NotAnId;
    }

    static std::pair<ID, bool> faceToEdge ( ID const& /*iFace*/, ID const& /*jEdge*/ )
    {
        ERROR_MSG ( "FaceToEdge not implemented for geo Line elements." );
        return std::make_pair ( static_cast<ID> (0), true );
    }
};


//! @ingroup BasRefSha
class Triangle
{
public:
    static const ReferenceShapes S_shape = TRIANGLE; //!< Identify the shape
    static const ReferenceGeometry S_geometry = FACE;//!< Identify the geometric entity
    static const UInt S_nDimensions = 2;             //!< Dimensionality
    static const UInt S_numVertices = 3;             //!< Number of vertices.
    static const UInt S_numEdges = 3;                //!< Number of edges
    static const UInt S_numFaces = 1;                //!< Number of faces
    static const UInt S_numFacets = S_numEdges;      //!< Number of facets
    static const UInt S_numRidges = S_numVertices;   //!< Number of ridges
    static const UInt S_numPeaks = 0;                //!< Number of peaks

    //! @return the local ID of the j-th point of the i-th face
    static ID faceToPoint ( ID const& /*iFace*/, ID const& jPoint )
    {
        return jPoint;
    };

    /*!
        @return a pair: the local numbering of the j-th edge on the i-th face, and
                true if the orientation of the edge on the face is consistent
                with that of the same edge on the element
    */
    static std::pair<ID, bool> faceToEdge ( ID const& /*iFace*/, ID const& jEdge )
    {
        return std::make_pair ( jEdge, true );
    }

    static ID facetToPeak ( ID const& /*iFacet*/, ID const& /*jPeak*/ )
    {
        return NotAnId;
    }
};


//! @ingroup BasRefSha
class Quad
{
public:
    static const ReferenceShapes S_shape = QUAD;     //!< Identify the shape
    static const ReferenceGeometry S_geometry = FACE;//!< Identify the geometric entity
    static const UInt S_nDimensions = 2;             //!< Dimensionality
    static const UInt S_numFaces = 1;                //!< Number of faces
    static const UInt S_numVertices = 4;             //!< Number of vertices.
    static const UInt S_numEdges = 4;                //!< Number of edges
    static const UInt S_numFacets = S_numEdges;      //!< Number of facets
    static const UInt S_numRidges = S_numVertices;   //!< Number of ridges
    static const UInt S_numPeaks = 0;                //!< Number of peaks

    //! @return the local ID of the j-th point of the i-th face
    static ID faceToPoint ( ID const& /*iFace*/, ID const& jPoint )
    {
        return jPoint;
    };

    /*!
        @return a pair: the local numbering of the j-th edge on the i-th face, and
                true if the orientation of the edge on the face is consistent
                with that of the same edge on the element
    */
    static std::pair<ID, bool> faceToEdge ( ID const& /*iFace*/, ID const& jEdge )
    {
        return std::make_pair ( jEdge, true );
    }

    static ID facetToPeak ( ID const& /*iFacet*/, ID const& /*jPeak*/ )
    {
        return NotAnId;
    }
};


//! @ingroup BasRefSha
class Tetra
{
public:
    static const ReferenceShapes S_shape = TETRA;             //!< Identify the shape
    static const ReferenceGeometry S_geometry = VOLUME;       //!< Identify the geometric entity
    static const UInt S_nDimensions = 3;                      //!< Dimensionality
    static const UInt S_numVertices = 4;                      //!< Number of vertices.
    static const UInt S_numFaces = 4;                         //!< Number of faces
    static const UInt S_numEdges = S_numFaces + S_numVertices - 2;//!< Number of edges
    static const UInt S_numFacets = S_numFaces;              //!< Number of facets
    static const UInt S_numRidges = S_numEdges;              //!< Number of ridges
    static const UInt S_numPeaks = S_numVertices;                        //!< Number of peaks

    /*!
        @return a pair: the local numbering of the j-th edge on the i-th face, and
                true if the orientation of the edge on the face is consistent
                with that of the same edge on the element
     */
    static std::pair<ID, bool> faceToEdge ( ID const& iFace, ID const& jEdge );
};


//! @ingroup BasRefSha
class Hexa
{
public:
    static const ReferenceShapes S_shape = HEXA;              //!< Identify the shape
    static const ReferenceGeometry S_geometry = VOLUME;       //!< Identify the geometric entity
    static const UInt S_nDimensions = 3;                      //!< Dimensionality
    static const UInt S_numFaces = 6;                         //!< Number of faces
    static const UInt S_numVertices = 8;                      //!< Number of vertices.
    static const UInt S_numEdges = S_numFaces + S_numVertices - 2;//!< Number of edges
    static const UInt S_numFacets = S_numFaces;              //!< Number of facets
    static const UInt S_numRidges = S_numEdges;              //!< Number of ridges
    static const UInt S_numPeaks = S_numVertices;                        //!< Number of peaks

    /*!
        @return a pair: the local numbering of the j-th edge on the i-th face, and
                true if the orientation of the edge on the face is consistent
                with that of the same edge on the element
     */
    static std::pair<ID, bool> faceToEdge ( ID const& iFace, ID const& jEdge );

};


/*
 Now the Basis Geometric Shapes, derived from the Reference Shapes
 */
//! @defgroup GeoShape Basis Geometric Shapes

// Forward Declarations.
class GeoPoint;
class LinearLine;
class QuadraticLine;
class LinearTriangle;
class QuadraticTriangle;
class LinearQuad;
class QuadraticQuad;
class LinearTetra;
class LinearTetraBubble;
class QuadraticTetra;
class LinearHexa;
class QuadraticHexa;


//! @ingroup GeoShape
//! A Geometric Shape
class GeoPoint:
    public Point
{
public:
    //! @name Public Types
    //@{
    typedef Point BasRefSha;
    typedef nullShape GeoBShape;             //!< Geometric shape of the boundary
    //@}
    static const UInt S_numPoints = 1; //!< Number of points
    static const UInt S_numPointsPerElement = 1;   //!< Number of points per element
    static const UInt S_numPointsPerFacet = 0;   //!< Number of points per facet
    static const UInt S_numPointsPerRidge = 0;   //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = 0;   //!< Number of points per peak

};


//! @ingroup GeoShape
//! A Geometric Shape
class LinearLine:
    public Line
{
public:
    //! @name Public Types
    //@{
    typedef Line BasRefSha;
    typedef GeoPoint GeoBShape;             //!< Geometric shape of the boundary
    //@}
    static const UInt S_numPoints = 2;      //!< Number of points
    static const UInt S_numPointsPerEdge = 0;   //!< Number of points per edge
    static const UInt S_numPointsPerVertex = 1; //!< Number of points per vertex
    static const UInt S_numPointsPerElement = S_numPointsPerEdge;   //!< Number of points per element
    static const UInt S_numPointsPerFacet = S_numPointsPerVertex;   //!< Number of points per facet
    static const UInt S_numPointsPerRidge = 0;   //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = 0;   //!< Number of points per peak
};
//! Inverts a line
template <>
inline ID reversePoint<LinearLine> ( ID const& pointId )
{
    static ID _rid[] = {1, 0};
    return _rid[pointId];
}


//! @ingroup GeoShape
//! A Geometric Shape
class QuadraticLine:
    public Line
{
public:
    //! @name Public Types
    //@{
    typedef Line BasRefSha;
    typedef GeoPoint GeoBShape;             //!< Geometric shape of the boundary
    //@}
    static const UInt S_numPoints = 3;      //!< Number of points
    static const UInt S_numPointsPerEdge = 1;   //!< Number of points per edge
    static const UInt S_numPointsPerVertex = 1; //!< Number of points per vertex
    static const UInt S_numPointsPerElement = S_numPointsPerEdge;   //!< Number of points per element
    static const UInt S_numPointsPerFacet = S_numPointsPerVertex;   //!< Number of points per facet
    static const UInt S_numPointsPerRidge = 0;   //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = 0;   //!< Number of points per peak

};

template <>
inline ID reversePoint<QuadraticLine> ( ID const& pointId )
{
    static ID _rid[] = {1, 0, 2};
    return _rid[pointId];
}


//! @ingroup GeoShape
//! A Geometric Shape
class LinearTriangle:
    public Triangle
{
public:
    //! @name Public Types
    //@{
    typedef Triangle BasRefSha;
    typedef LinearLine GeoBShape;
    //@}
    static const UInt S_numPoints = 3;      //!< Number of points
    static const UInt S_numPointsPerVertex = 1; //!< Number of points per vertex
    static const UInt S_numPointsPerEdge = 0;   //!< Number of points per edge
    static const UInt S_numPointsPerFace = 0;   //!< Number of points per face
    static const UInt S_numPointsPerElement = S_numPointsPerFace;  //!< Number of points per element
    static const UInt S_numPointsPerFacet = S_numPointsPerEdge;  //!< Number of points per facet
    static const UInt S_numPointsPerRidge = S_numPointsPerVertex;   //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = 0;   //!< Number of points per peak

    //! @return the local ID of the j-th point of the i-th edge
    static ID edgeToPoint ( ID const& iEdge, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th facet
    static ID facetToPoint ( ID const& iFacet, ID const& jPoint )
    {
        return edgeToPoint ( iFacet, jPoint );
    }

    //! @return the local ID of the j-th ridge of the i-th facet
    static ID facetToRidge ( ID const& iFacet, ID const& jRidge )
    {
        return edgeToPoint (iFacet, jRidge);
    }
};

template <>
inline ID reversePoint<LinearTriangle> ( ID const& pointId )
{
    static ID _rid[] = {1, 0, 2};
    return _rid[pointId];
}


//! @ingroup GeoShape
//! A Geometric Shape
class QuadraticTriangle:
    public Triangle
{
public:
    //! @name Public Types
    //@{
    typedef Triangle BasRefSha;
    typedef QuadraticLine GeoBShape;
    //@}
    static const UInt S_numPoints = 6;      //!< Number of points
    static const UInt S_numPointsPerVertex = 1; //!< Number of points per vertex
    static const UInt S_numPointsPerEdge = 1;   //!< Number of points per edge
    static const UInt S_numPointsPerFace = 0;   //!< Number of points per face
    static const UInt S_numPointsPerElement = S_numPointsPerFace;  //!< Number of points per element
    static const UInt S_numPointsPerFacet = S_numPointsPerEdge;  //!< Number of points per facet
    static const UInt S_numPointsPerRidge = S_numPointsPerVertex; //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = 0;   //!< Number of points per peak


    //! @return the local ID of the j-th point of the i-th edge
    static ID edgeToPoint ( ID const& iEdge, ID const& jPoint );
    //! @return the local ID of the j-th point of the i-th facet
    static ID facetToPoint ( ID const& iFacet, ID const& jPoint )
    {
        return edgeToPoint ( iFacet, jPoint );
    }
    static ID facetToRidge ( ID const& iFacet, ID const& jRidge )
    {
        return edgeToPoint (iFacet, jRidge);
    }
};

template <>
inline ID reversePoint<QuadraticTriangle> ( ID const& pointId )
{
    static ID _rid[] = {1, 0, 2, 3, 5, 4};
    return _rid[pointId];
}


//! @ingroup GeoShape
//! A Geometric Shape
class LinearQuad:
    public Quad
{
public:
    //! @name Public Types
    //@{
    typedef Quad BasRefSha;
    typedef LinearLine GeoBShape;
    //@}
    static const UInt S_numPoints = 4;      //!< Number of points
    static const UInt S_numPointsPerVertex = 1; //!< Number of points per vertex
    static const UInt S_numPointsPerEdge = 0;   //!< Number of points per edge
    static const UInt S_numPointsPerFace = 0;   //!< Number of points per face
    static const UInt S_numPointsPerElement = S_numPointsPerFace;  //!< Number of points per element
    static const UInt S_numPointsPerFacet = S_numPointsPerEdge;  //!< Number of points per facet
    static const UInt S_numPointsPerRidge = S_numPointsPerVertex;  //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = 0;   //!< Number of points per peak


    //! @return the local ID of the j-th point of the i-th edge
    static ID edgeToPoint ( ID const& iEdge, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th facet
    static ID facetToPoint ( ID const& iFacet, ID const& jPoint )
    {
        return edgeToPoint ( iFacet, jPoint );
    }

    //! @return the local ID of the j-th point of the i-th ridge
    static ID facetToRidge ( ID const& iFacet, ID const& jRidge )
    {
        return edgeToPoint (iFacet, jRidge);
    }
};

//! Specialization
template <>
inline ID reversePoint<LinearQuad> ( ID const& pointId )
{
    static ID _rid[] = {3, 2, 1, 0};
    return _rid[pointId];
}


//! @ingroup GeoShape
//! A Geometric Shape
class QuadraticQuad:
    public Quad
{
public:
    //! @name Public Types
    //@{
    typedef Quad BasRefSha;
    typedef QuadraticLine GeoBShape;
    //@}
    static const UInt S_numPoints = 9;      //!< Number of points
    static const UInt S_numPointsPerVertex = 1; //!< Number of points per vertex
    static const UInt S_numPointsPerEdge = 1;   //!< Number of points per edge
    static const UInt S_numPointsPerFace = 1;   //!< Number of points per face
    static const UInt S_numPointsPerElement = S_numPointsPerFace;  //!< Number of points per element
    static const UInt S_numPointsPerFacet = S_numPointsPerEdge;  //!< Number of points per facet
    static const UInt S_numPointsPerRidge = S_numPointsPerVertex;   //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = 0;   //!< Number of points per peak

    //! @return the local ID of the j-th point of the i-th edge
    static ID edgeToPoint ( ID const& iEdge, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th facet
    static ID facetToPoint ( ID const& iFacet, ID const& jPoint )
    {
        return edgeToPoint ( iFacet, jPoint );
    }

    //! @return the local ID of the j-th point of the i-th ridge
    static ID facetToRidge ( ID const& iFacet, ID const& jRidge )
    {
        return edgeToPoint (iFacet, jRidge);
    }
};
//! Specialization
template <>
inline ID reversePoint<QuadraticQuad> ( ID const& pointId )
{
    static ID _rid[] = {3, 2, 1, 0, 6, 5, 4, 7, 8};
    return _rid[pointId];
}


//! @ingroup GeoShape
//! A Geometric Shape
class LinearTetra:
    public Tetra
{
public:
    //! @name Public Types
    //@{
    typedef Tetra BasRefSha;
    typedef LinearTriangle GeoBShape;
    //@}
    static const UInt S_numPoints = 4;      //!< Number of points
    static const UInt S_numPointsPerVertex = 1; //!< Number of points per vertex
    static const UInt S_numPointsPerEdge = 0;   //!< Number of points per edge
    static const UInt S_numPointsPerFace = 0;   //!< Number of points per face
    static const UInt S_numPointsPerVolume = 0; //!< Number of points per volume
    static const UInt S_numPointsPerElement = S_numPointsPerVolume; //!< Number of points per element
    static const UInt S_numPointsPerFacet =  S_numPointsPerFace; //!< Number of points per facet
    static const UInt S_numPointsPerRidge = S_numPointsPerEdge;    //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = S_numPointsPerVertex;   //!< Number of points per peak

    //! @return the local ID of the j-th point of the i-th edge
    static ID edgeToPoint ( ID const& iEdge, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th face
    static ID faceToPoint ( ID const& iFace, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th facet
    static ID facetToPoint ( ID const& iFacet, ID const& jPoint )
    {
        return faceToPoint ( iFacet, jPoint );
    }

    //! @return the local ID of the j-th ridge of the i-th facet
    static ID facetToRidge ( ID const& iFacet, ID const& jRidge )
    {
        return faceToEdge (iFacet, jRidge).first;
    }

    //! @return the local ID of the j-th peak of the i-th facet
    static ID facetToPeak ( ID const& iFacet, ID const& jPeak )
    {
        return faceToPoint (iFacet, jPeak);
    }
};

//! Specialization
template <>
inline ID reversePoint<LinearTetra> ( ID const& pointId )
{
    static ID _rid[] = {1, 0, 2, 3};
    return _rid[pointId];
}


//! @ingroup GeoShape
//! A Geometric Shape
class LinearTetraBubble:
    public Tetra
{
public:
    //! @name Public Types
    //@{
    typedef Tetra BasRefSha;
    typedef LinearTriangle GeoBShape;
    //@}
    static const UInt S_numPoints = 5;      //!< Number of points
    static const UInt S_numPointsPerVertex = 1; //!< Number of points per vertex
    static const UInt S_numPointsPerEdge = 0;   //!< Number of points per edge
    static const UInt S_numPointsPerFace = 0;   //!< Number of points per face
    static const UInt S_numPointsPerVolume = 1; //!< Number of points per volume
    static const UInt S_numPointsPerElement = S_numPointsPerVolume; //!< Number of points per element
    static const UInt S_numPointsPerFacet =  S_numPointsPerFace; //!< Number of points per facet
    static const UInt S_numPointsPerRidge = S_numPointsPerEdge;    //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = S_numPointsPerVertex;   //!< Number of points per peak

    //! @return the local ID of the j-th point of the i-th edge
    static ID edgeToPoint ( ID const& iEdge, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th face
    static ID faceToPoint ( ID const& iFace, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th facet
    static ID facetToPoint ( ID const& iFacet, ID const& jPoint )
    {
        return faceToPoint ( iFacet, jPoint );
    }

    //! @return the local ID of the j-th ridge of the i-th facet
    static ID facetToRidge ( ID const& iFacet, ID const& jRidge )
    {
        return faceToEdge (iFacet, jRidge).first;
    }

    //! @return the local ID of the j-th peak of the i-th facet
    static ID facetToPeak ( ID const& iFacet, ID const& jPeak )
    {
        return faceToPoint (iFacet, jPeak);
    }
};

template <>
inline ID reversePoint<LinearTetraBubble> ( ID const& pointId )
{
    static ID _rid[] = {1, 0, 2, 3, 4};
    return _rid[pointId];
}


//! @ingroup GeoShape
//! A Geometric Shape
class QuadraticTetra:
    public Tetra
{
public:
    //! @name Public Types
    //@{
    typedef Tetra BasRefSha;
    typedef QuadraticTriangle GeoBShape;
    //@}
    static const UInt S_numPoints = 10;     //!< Number of points
    static const UInt S_numPointsPerVertex = 1; //!< Number of points per vertex
    static const UInt S_numPointsPerEdge = 1;   //!< Number of points per edge
    static const UInt S_numPointsPerFace = 0;   //!< Number of points per face
    static const UInt S_numPointsPerVolume = 0; //!< Number of points per volume
    static const UInt S_numPointsPerElement = S_numPointsPerVolume; //!< Number of points per element
    static const UInt S_numPointsPerFacet =  S_numPointsPerFace; //!< Number of points per facet
    static const UInt S_numPointsPerRidge = S_numPointsPerEdge;    //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = S_numPointsPerVertex;   //!< Number of points per peak


    //! @return the local ID of the j-th point of the i-th edge
    static ID edgeToPoint ( ID const& iEdge, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th face
    static ID faceToPoint ( ID const& iFace, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th facet
    static ID facetToPoint ( ID const& iFacet, ID const& jPoint )
    {
        return faceToPoint ( iFacet, jPoint );
    }

    //! @return the local ID of the j-th ridge of the i-th facet
    static ID facetToRidge ( ID const& iFacet, ID const& jRidge )
    {
        return faceToEdge (iFacet, jRidge).first;
    }

    //! @return the local ID of the j-th peak of the i-th facet
    static ID facetToPeak ( ID const& iFacet, ID const& jPeak )
    {
        return faceToPoint (iFacet, jPeak);
    }

};

template <>
inline ID reversePoint<QuadraticTetra> ( ID const& pointId )
{
    static ID _rid[] = {1, 0, 2, 3, 4, 6, 5, 8, 7, 9};
    return _rid[pointId];
}


//! @ingroup GeoShape
//! A Geometric Shape
class LinearHexa:
    public Hexa
{
public:
    //! @name Public Types
    //@{
    typedef Hexa BasRefSha;
    typedef LinearQuad GeoBShape;
    //@}
    static const UInt S_numPoints = 8;      //!< Number of points
    static const UInt S_numPointsPerVertex = 1; //!< Number of points per vertex
    static const UInt S_numPointsPerEdge = 0;   //!< Number of points per edge
    static const UInt S_numPointsPerFace = 0;   //!< Number of points per face
    static const UInt S_numPointsPerVolume = 0; //!< Number of points per volume
    static const UInt S_numPointsPerElement = S_numPointsPerVolume; //!< Number of points per element
    static const UInt S_numPointsPerFacet =  S_numPointsPerFace; //!< Number of points per facet
    static const UInt S_numPointsPerRidge = S_numPointsPerEdge;    //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = S_numPointsPerVertex;   //!< Number of points per peak


    //! @return the local ID of the j-th point of the i-th edge
    static ID edgeToPoint ( ID const& iEdge, ID const& jPoint );
    //! @return the local ID of the j-th point of the i-th face
    static ID faceToPoint ( ID const& iFace, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th facet
    static ID facetToPoint ( ID const& iFacet, ID const& jPoint )
    {
        return faceToPoint ( iFacet, jPoint );
    }

    //! @return the local ID of the j-th ridge of the i-th facet
    static ID facetToRidge ( ID const& iFacet, ID const& jRidge )
    {
        return faceToEdge (iFacet, jRidge).first;
    }

    //! @return the local ID of the j-th peak of the i-th facet
    static ID facetToPeak ( ID const& iFacet, ID const& jPeak )
    {
        return faceToPoint (iFacet, jPeak);
    }

};

template <>
inline ID reversePoint<LinearHexa> ( ID const& pointId )
{
    static ID _rid[] = {3, 2, 1, 0, 7, 6, 5, 4};
    return _rid[pointId];
}


//! @ingroup GeoShape
//! A Geometric Shape
class QuadraticHexa:
    public Hexa
{
public:
    //! @name Public Types
    //@{
    typedef Hexa BasRefSha;
    typedef QuadraticQuad GeoBShape;
    //@}
    static const UInt S_numPoints = 27;     //!< Number of points
    static const UInt S_numPointsPerVertex = 1; //!< Number of points per vertex
    static const UInt S_numPointsPerEdge = 1;   //!< Number of points per edge
    static const UInt S_numPointsPerFace = 1;   //!< Number of points per face
    static const UInt S_numPointsPerVolume = 1; //!< Number of points per volume
    static const UInt S_numPointsPerElement = S_numPointsPerVolume; //!< Number of points per element
    static const UInt S_numPointsPerFacet =  S_numPointsPerFace; //!< Number of points per facet
    static const UInt S_numPointsPerRidge = S_numPointsPerEdge;    //!< Number of points per ridge
    static const UInt S_numPointsPerPeak = S_numPointsPerVertex;   //!< Number of points per peak


    //! @return the local ID of the j-th point of the i-th edge
    static ID edgeToPoint ( ID const& iEdge, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th face
    static ID faceToPoint ( ID const& iFace, ID const& jPoint );

    //! @return the local ID of the j-th point of the i-th facet
    static ID facetToPoint ( ID const& iFacet, ID const& jPoint )
    {
        return faceToPoint ( iFacet, jPoint );
    }

    //! @return the local ID of the j-th ridge of the i-th facet
    static ID facetToRidge ( ID const& iFacet, ID const& jRidge )
    {
        return faceToEdge (iFacet, jRidge).first;
    }

    //! @return the local ID of the j-th peak of the i-th facet
    static ID facetToPeak ( ID const& iFacet, ID const& jPeak )
    {
        return faceToPoint (iFacet, jPeak);
    }

};

template <>
inline ID reversePoint<QuadraticHexa> ( ID const& pointId )
{
    static ID _rid[] = {3, 2, 1, 0, 7, 6, 5, 4,
                        10, 9, 8, 11, 15, 14, 13, 12,
                        18, 17, 16, 19, 20, 23, 22, 21, 24, 25, 26
                       };
    return _rid[pointId];
}


}

#endif  // ELEMENTSHAPES_H

