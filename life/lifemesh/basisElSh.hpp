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
	<li> <tt>static ID eToP(ID const jEdge, ID const iPoint);</tt>
         which returns the local ID of the i-th point of the j-th edge of the GeoShape

	<li> <tt>static ID fToP(ID const jFace, ID const iPoint);</tt>
         which returns the local ID of the i-th point of the j-th face of the GeoShape

	<li> <tt>static pair<ID,bool> fToE(ID const jFace, ID const iEdge);</tt>
	     which returns the local numbering of the i-th edge on the j-th face.
         It returns also if the orientation of the edge on the face is consistent
	     with that of the same edge on the element
	</ol>

	@note The methods edge to point-ID (EtoP) and Face to point-ID (FtoP)
	      return the local id number of points on faces and edges (when relevant)
	@note We follow the convention of indicating THE VERTICES FIRST in the list
	      of dofs
*/

#ifndef BASISELSH_H
#define BASISELSH_H 1

#include <utility>
#include <life/lifecore/life.hpp>

namespace LifeV
{
//! A utility to invert point numbering on a GeoShape
template <typename GeoShapeType>
class reversePoint
{
public:
    //! @name Operators
    //@{

    //! The function call operator
	/*!
	    <ol>
	    <li> If the point Id is smaller than GeoShapeType::numVertices, return<br>
    		GeoShapeType::numVertices - pointId + 1
	    <li> If the point Id is bigger than GeoShapeType::numVertices, return<br>
            GeoShapeType::numPoints - point + GeoShapeType::numVertices + 1;
	    </ol>

        @note this numbering follows the convention that VERTICES are always numbered first.
        @param pointId the old point ID
        @return the new point ID
	 */
    inline ID operate( ID const & pointId );
    //@}
};


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
UInt getReferenceShapeDimension(const ReferenceShapes& shape);

UInt __attribute__ ((__deprecated__)) getReferenceDimension(const ReferenceShapes& shape);


/*! @enum ReferenceShapes
    Lists of the geometries of the finite element in the reference space,
    supported by the library

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
class Point
{
public:
    static const ReferenceShapes Shape      = POINT; //!< Identify the shape
    static const ReferenceGeometry Geometry = VERTEX;//!< Identify the geometric entity
    static const UInt nDim                  = 0;     //!< Dimensionality
    static const UInt numFaces              = 0;     //!< Number of faces
    static const UInt numEdges              = 0;     //!< Number of edges
    static const UInt numVertices           = 1;     //!< Number of vertices.
};


//! @ingroup BasRefSha
class Line
{
public:
    static const ReferenceShapes Shape = LINE;     //!< Identify the shape
    static const ReferenceGeometry Geometry = EDGE;//!< Identify the geometric entity
    static const UInt nDim = 1;                    //!< Dimensionality
    static const UInt numFaces = 0;                //!< Number of faces
    static const UInt numEdges = 1;                //!< Number of edges
    static const UInt numVertices = 2;             //!< Number of vertices.
};


//! @ingroup BasRefSha
class Triangle
{
public:
    static const ReferenceShapes Shape = TRIANGLE; //!< Identify the shape
    static const ReferenceGeometry Geometry = FACE;//!< Identify the geometric entity
    static const UInt nDim = 2;                    //!< Dimensionality
    static const UInt numVertices = 3;             //!< Number of vertices.
    static const UInt numFaces = 1;                //!< Number of faces
    static const UInt numEdges = numVertices;      //!< Number of edges
};


//! @ingroup BasRefSha
class Quad
{
public:
    static const ReferenceShapes Shape = QUAD;     //!< Identify the shape
    static const ReferenceGeometry Geometry = FACE;//!< Identify the geometric entity
    static const UInt nDim = 2;                    //!< Dimensionality
    static const UInt numFaces = 1;                //!< Number of faces
    static const UInt numVertices = 4;             //!< Number of vertices.
    static const UInt numEdges = numVertices;      //!< Number of edges
};


//! @ingroup BasRefSha
class Tetra
{
public:
    static const ReferenceShapes Shape = TETRA;             //!< Identify the shape
    static const ReferenceGeometry Geometry = VOLUME;       //!< Identify the geometric entity
    static const UInt nDim = 3;                             //!< Dimensionality
    static const UInt numVertices = 4;                      //!< Number of vertices.
    static const UInt numFaces = 4;                         //!< Number of faces
    static const UInt numEdges = numFaces + numVertices - 2;//!< Number of edges
};


//! @ingroup BasRefSha
class Hexa
{
public:
    static const ReferenceShapes Shape = HEXA;              //!< Identify the shape
    static const ReferenceGeometry Geometry = VOLUME;       //!< Identify the geometric entity
    static const UInt nDim = 3;                             //!< Dimensionality
    static const UInt numFaces = 6;                         //!< Number of faces
    static const UInt numVertices = 8;                      //!< Number of vertices.
    static const UInt numEdges = numFaces + numVertices - 2;//!< Number of edges
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
    //@}
    static const UInt numPoints = 1; //!< Number of points
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
    typedef GeoPoint GeoBShape;           //!< Geometric shape of the boundary
    //@}
    static const UInt numPoints = 2;      //!< Number of points
    static const UInt nbPtsPerVertex = 1; //!< Number of points per vertex
    static const UInt nbPtsPerEdge = 0;   //!< Number of points per edge
};


//! @ingroup GeoShape
//! A Geometric Shape
class QuadraticLine:
        public Line
{
public:
    //! @name Public Types
    //@{
    typedef Line BasRefSha;
    typedef GeoPoint GeoBShape;           //!< Geometric shape of the boundary
    //@}
    static const UInt numPoints = 3;      //!< Number of points
    static const UInt nbPtsPerVertex = 1; //!< Number of points per vertex
    static const UInt nbPtsPerEdge = 1; //!< Number of points per edge
};


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
    static const UInt numPoints = 3; //!< Number of points
    static const UInt nbPtsPerVertex = 1; //!< Number of points per vertex
    static const UInt nbPtsPerEdge = 0; //!< Number of points per edge
    static const UInt nbPtsPerFace = 0; //!< Number of points per face
    //! @return the local ID of the j-th point of the i-th edge
    static ID eToP( ID const iEdge, ID const jPoint );
};


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
    static const UInt numPoints = 6; //!< Number of points
    static const UInt nbPtsPerVertex = 1; //!< Number of points per vertex
    static const UInt nbPtsPerEdge = 1; //!< Number of points per edge
    static const UInt nbPtsPerFace = 0; //!< Number of points per face
    //! @return the local ID of the j-th point of the i-th edge
    static ID eToP( ID const iEdge, ID const jPoint );
};


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
    static const UInt numPoints = 4; //!< Number of points
    static const UInt nbPtsPerVertex = 1; //!< Number of points per vertex
    static const UInt nbPtsPerEdge = 0; //!< Number of points per edge
    static const UInt nbPtsPerFace = 0; //!< Number of points per face
    //! @return the local ID of the j-th point of the i-th edge
    static ID eToP( ID const iEdge, ID const jPoint );
};


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
    static const UInt numPoints = 9; //!< Number of points
    static const UInt nbPtsPerVertex = 1; //!< Number of points per vertex
    static const UInt nbPtsPerEdge = 1; //!< Number of points per edge
    static const UInt nbPtsPerFace = 1; //!< Number of points per face
    //! @return the local ID of the j-th point of the i-th edge
    static ID eToP( ID const iEdge, ID const jPoint );
};


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
    static const UInt numPoints = 4; //!< Number of points
    static const UInt nbPtsPerVertex = 1; //!< Number of points per vertex
    static const UInt nbPtsPerEdge = 0; //!< Number of points per edge
    static const UInt nbPtsPerFace = 0; //!< Number of points per face
    static const UInt nbPtsPerVolume = 0; //!< Number of points per volume
    //! @return the local ID of the j-th point of the i-th edge
    static ID eToP( ID const iEdge, ID const jPoint );
    //! @return the local ID of the j-th point of the i-th face
    static ID fToP( ID const iFace, ID const jPoint );
    /*!
        @return a pair: the local numbering of the j-th edge on the i-th face, and
                true if the orientation of the edge on the face is consistent
	            with that of the same edge on the element
     */
    static std::pair<ID, bool> fToE( ID const iFace, ID const jEdge );
};


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
    static const UInt numPoints = 5; //!< Number of points
    static const UInt nbPtsPerVertex = 1; //!< Number of points per vertex
    static const UInt nbPtsPerEdge = 0; //!< Number of points per edge
    static const UInt nbPtsPerFace = 0; //!< Number of points per face
    static const UInt nbPtsPerVolume = 1; //!< Number of points per volume
    //! @return the local ID of the j-th point of the i-th edge
    static ID eToP( ID const iEdge, ID const jPoint );
    //! @return the local ID of the j-th point of the i-th face
    static ID fToP( ID const iFace, ID const jPoint );
    /*!
        @return a pair: the local numbering of the j-th edge on the i-th face, and
                true if the orientation of the edge on the face is consistent
	            with that of the same edge on the element
     */
    static std::pair<ID, bool> fToE( ID const iFace, ID const jEdge );
};


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
    static const UInt numPoints = 10; //!< Number of points
    static const UInt nbPtsPerVertex = 1; //!< Number of points per vertex
    static const UInt nbPtsPerEdge = 1; //!< Number of points per edge
    static const UInt nbPtsPerFace = 0; //!< Number of points per face
    static const UInt nbPtsPerVolume = 0; //!< Number of points per volume
    //! @return the local ID of the j-th point of the i-th edge
    static ID eToP( ID const iEdge, ID const jPoint );
    //! @return the local ID of the j-th point of the i-th face
    static ID fToP( ID const iFace, ID const jPoint );
    /*!
        @return a pair: the local numbering of the j-th edge on the i-th face, and
                true if the orientation of the edge on the face is consistent
	            with that of the same edge on the element
     */
    static std::pair<ID, bool> fToE( ID const iFace, ID const jEdge );
};


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
    static const UInt numPoints = 8; //!< Number of points
    static const UInt nbPtsPerVertex = 1; //!< Number of points per vertex
    static const UInt nbPtsPerEdge = 0; //!< Number of points per edge
    static const UInt nbPtsPerFace = 0; //!< Number of points per face
    static const UInt nbPtsPerVolume = 0; //!< Number of points per volume
    //! @return the local ID of the j-th point of the i-th edge
    static ID eToP( ID const iEdge, ID const jPoint );
    //! @return the local ID of the j-th point of the i-th face
    static ID fToP( ID const iFace, ID const jPoint );
    /*!
        @return a pair: the local numbering of the j-th edge on the i-th face, and
                true if the orientation of the edge on the face is consistent
	            with that of the same edge on the element
     */
    static std::pair<ID, bool> fToE( ID const iFace, ID const jEdge );
};


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
    static const UInt numPoints = 27; //!< Number of points
    static const UInt nbPtsPerVertex = 1; //!< Number of points per vertex
    static const UInt nbPtsPerEdge = 1; //!< Number of points per edge
    static const UInt nbPtsPerFace = 1; //!< Number of points per face
    static const UInt nbPtsPerVolume = 1; //!< Number of points per volume
    //! @return the local ID of the j-th point of the i-th edge
    static ID eToP( ID const iEdge, ID const jPoint );
    //! @return the local ID of the j-th point of the i-th face
    static ID fToP( ID const iFace, ID const jPoint );
    /*!
        @return a pair: the local numbering of the j-th edge on the i-th face, and
                true if the orientation of the edge on the face is consistent
	            with that of the same edge on the element
     */
    static std::pair<ID, bool> fToE( ID const iFace, ID const jEdge );
};


/*******************************************************************
          IMPLEMENTATION
*******************************************************************/

// ===================================================
// Operators
// ===================================================
template <typename GeoShapeType>
inline
ID reversePoint<GeoShapeType>::
operate( ID const & point )
{
    return point <= GeoShapeType::numVertices ?
    		GeoShapeType::numVertices - point + 1 :
            GeoShapeType::numPoints - point + GeoShapeType::numVertices + 1;
}
}

#endif

