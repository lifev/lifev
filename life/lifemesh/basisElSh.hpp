//! \file basisElSh.h
/*! Contains the basic element shapes, to be used by Geometric and Finite
  Elements 

  $Header: /cvsroot/lifev/lifev/life/lifemesh/basisElSh.hpp,v 1.2 2004-02-24 13:28:11 prudhomm Exp $

 \version 0.0 Experimental   19/8/99. Luca Formaggia

 These classes provides the most basic definitions for the geometric
 elements. The variables names are self-explanatory.

 The Basis Geometric Shapes (<tt>GeoShape</tt>) are derived from a set of
 class called Basis Reference Shapes (<tt>BasRefSha</tt>), which contains
 the very basic information about the reference element geometries
 supported by the library (reference element geometry=geometry of the
 finite element in the reference space.)

 A Basis Geometric Shape contains (when relevant) the following methods

  static ID eToP(ID const _localEdge, ID const _point);


  static ID fToP(ID const _localFace, ID const _point);

  which returns the local ID of a point on a face of the GeoShape

  static pair<ID,bool> fToE(ID const _localFace, ID const _edge);

  which returns the local ID corresponding and edge on the face _localface.
  It returna also if the orientation of the edge on the face is consistent
  with that of the same edge on the element
  

 \note The methods edge to point-ID (EtoP) and Face to point-ID (FtoP)
 return the local id number of points on faces and edges (when relevant)
 We follow the convention of indicating THE VERTICES FIRST!  */


#ifndef _BASISELSH_HH_
#define _BASISELSH_HH_
#if defined(__Linux)
//# include <pair.h>
# include <utility>
#elif defined(__OSF1)
# include <utility>
#endif

#include "lifeV.hpp"
using std::pair;
using std::make_pair;


//! An utility to  invert point numbering on a GeoShape
template <typename GeoShape>
class reversePoint
{
 public:
  INLINE ID operate(ID const & point);
};


//                   *********** BASIS REFERENCE SHAPES ****************

enum ReferenceShapes {NONE, POINT, LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA};
enum ReferenceGeometry {VERTEX=0,EDGE=1,FACE=2,VOLUME=3};
//! \defgroup BasRefSha   Basis Reference Shapes

//! \ingroup BasRefSha
class Point
{
public:
  static const ReferenceShapes Shape=POINT; //!< Identify the shape
  static const ReferenceGeometry Geometry=VERTEX;
  static const UInt nDim=0;   //!< Dimensionality  
  static const UInt numFaces=0;//!< Number of faces
  static const UInt numEdges=0;//!< Number of faces
  static const UInt numVertices=1; //!< Number of vertices.
};

//! \ingroup BasRefSha
class Line
{
 public:
  static const ReferenceShapes Shape=LINE;
  static const ReferenceGeometry Geometry=EDGE;
  static const UInt nDim=1;
  static const UInt numFaces=0;
  static const UInt numEdges=0;
  static const UInt numVertices=2;
};

//! \ingroup BasRefSha
class Triangle
{
 public:
  static const ReferenceShapes Shape=TRIANGLE;
  static const ReferenceGeometry Geometry=FACE;
  static const UInt nDim=2;
  static const UInt numVertices=3;
  static const UInt numFaces=0;
  static const UInt numEdges=numVertices;
};

//! \ingroup BasRefSha
class Quad
{
 public:
  static const ReferenceShapes Shape=QUAD;
  static const ReferenceGeometry Geometry=FACE;
  static const UInt nDim=2;
  static const UInt numFaces=0;
  static const UInt numVertices=4;
  static const UInt numEdges=numVertices;
};

//! \ingroup BasRefSha
class Tetra
{
 public:
  static const ReferenceShapes Shape=TETRA;
  static const ReferenceGeometry Geometry=VOLUME;
  static const UInt nDim=3;
  static const UInt numVertices=4;
  static const UInt numFaces=4;
  static  const UInt numEdges=numFaces+numVertices-2;
};

//! \ingroup BasRefSha
class Hexa
{
 public:
  static const ReferenceShapes Shape=HEXA;
  static const ReferenceGeometry Geometry=VOLUME;
  static const UInt nDim=3;
  static const UInt numFaces=6;
  static const UInt numVertices=8;
  static  const UInt numEdges=numFaces+numVertices-2;
};
/*
 Now the Basis Geometric Shapes, derived from the Reference Shapes
*/
//! \defgroup GeoShape Basis Geometric Shapes
/// Forward Declarations.
class GeoPoint;
class LinearLine;
class QuadraticLine;
class LinearTriangle;
class QuadraticTriangle;
class LinearQuad;
class QuadraticQuad;
class LinearTetra;
class QuadraticTetra;
class LinearHexa;
class QuadraticHexa;

//! \ingroup GeoShape
class GeoPoint: 
public Point
{
public: 
  typedef Point BasRefSha;
  static const UInt numPoints=1;//!< Number of points
};

//! \ingroup GeoShape
class LinearLine:
public Line
{
public:
  typedef Line BasRefSha;
  typedef GeoPoint GeoBShape;//!< GeoShape of the boundary
  static const UInt numPoints=2;//!< Number of points
  static const UInt  nbPtsPerVertex  = 1;//!< Number of points per vertex in this GeoShape
  static const UInt  nbPtsPerEdge    = 0;//!< Number of points per edge in this GeoShape
};

//! \ingroup GeoShape
class QuadraticLine:
public Line
{
public:
  typedef Line BasRefSha;
  typedef GeoPoint GeoBShape;
  static const UInt numPoints=3;
  static const UInt  nbPtsPerVertex  = 1;
  static const UInt  nbPtsPerEdge    = 1;
};


//! \ingroup GeoShape
class LinearTriangle:
public Triangle
{
public:
  typedef Triangle BasRefSha;
  typedef LinearLine GeoBShape;
  static const UInt numPoints=3;
  static const UInt  nbPtsPerVertex  = 1;
  static const UInt  nbPtsPerEdge    = 0;
  static const UInt  nbPtsPerFace    = 0;
  static ID eToP(ID const _localEdge, ID const _point); 
};

//! \ingroup GeoShape
class QuadraticTriangle:
public Triangle
{
public:
  typedef Triangle BasRefSha;
  typedef QuadraticLine GeoBShape;
  static const UInt numPoints=6;
  static const UInt  nbPtsPerVertex  = 1;
  static const UInt  nbPtsPerEdge    = 1;
  static const UInt  nbPtsPerFace    = 0;
  static ID eToP(ID const _localEdge, ID const  _point);
};

//! \ingroup GeoShape
class LinearQuad:
public Quad
{
public:
  typedef Quad BasRefSha;
  typedef LinearLine GeoBShape;
  static const UInt numPoints=4;
  static const UInt  nbPtsPerVertex  = 1;
  static const UInt  nbPtsPerEdge    = 0;
  static const UInt  nbPtsPerFace    = 0;
  static ID eToP(ID const _localEdge, ID const  _point);
};

//! \ingroup GeoShape
class QuadraticQuad:
public Quad
{
public:
  typedef Quad BasRefSha;
  typedef QuadraticLine GeoBShape;
  static const UInt numPoints=9;
  static const UInt  nbPtsPerVertex  = 1;
  static const UInt  nbPtsPerEdge    = 1;
  static const UInt  nbPtsPerFace    = 1;
  static ID eToP(ID const _localEdge, ID const  _point);
};

//! \ingroup GeoShape
class LinearTetra:
public Tetra
{
public:
  typedef Tetra BasRefSha;
  typedef LinearTriangle GeoBShape;
  static const UInt numPoints=4;
  static const UInt  nbPtsPerVertex  = 1;
  static const UInt  nbPtsPerEdge    = 0;
  static const UInt  nbPtsPerFace    = 0;
  static const UInt  nbPtsPerVolume  = 0;
  static ID eToP(ID const _localEdge, ID const _point);
  static ID fToP(ID const _localFace, ID const _point);
  static pair<ID,bool> fToE(ID const _localFace, ID const _edge);
};

//! \ingroup GeoShape
class QuadraticTetra:
public Tetra
{
public:
  typedef Tetra BasRefSha;
  typedef QuadraticTriangle GeoBShape;
  static const UInt numPoints=10;
  static const UInt  nbPtsPerVertex  = 1;
  static const UInt  nbPtsPerEdge    = 1;
  static const UInt  nbPtsPerFace    = 0;
  static const UInt  nbPtsPerVolume  = 0;
  static ID eToP(ID const _localEdge, ID const _point);
  static ID fToP(ID const _localFace, ID const _point);
  static pair<ID,bool> fToE(ID const _localFace, ID const _edge);
};

//! \ingroup GeoShape
class LinearHexa:
public Hexa
{
public:
  typedef Hexa BasRefSha;
  typedef LinearQuad GeoBShape;
  static const UInt numPoints=8;
  static const UInt  nbPtsPerVertex  = 1;
  static const UInt  nbPtsPerEdge    = 0;
  static const UInt  nbPtsPerFace    = 0;
  static const UInt  nbPtsPerVolume  = 0;
  static ID eToP(ID const _localEdge, ID const _point);
  static ID fToP(ID const _localFace, ID const _point);
  static pair<ID,bool> fToE(ID const _localFace, ID const _edge);
};

//! \ingroup GeoShape
class QuadraticHexa:
public Hexa
{
public:
  typedef Hexa BasRefSha;
  typedef QuadraticQuad GeoBShape;
  static const UInt numPoints=27;
  static const UInt  nbPtsPerVertex  = 1;
  static const UInt  nbPtsPerEdge    = 1;
  static const UInt  nbPtsPerFace    = 1;
  static const UInt  nbPtsPerVolume  = 1;
  static ID eToP(ID const _localEdge, ID const _point);
  static ID fToP(ID const _localFace, ID const _point);
  static pair<ID,bool> fToE(ID const _localFace, ID const _edge);
};

/*******************************************************************
          IMPLEMENTATION
*******************************************************************/
template <typename GeoShape>
INLINE
ID reversePoint<GeoShape>::
operate(ID const & point)
{
  return point <=GeoShape::numVertices? GeoShape::numVertices-point+1:GeoShape::numPoints-point+GeoShape::numVertices+1;
};


#endif

