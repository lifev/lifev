/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/*!
  \file basicOneDMesh.hpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief This file contains a very basic one d mesh handler.

*/

#ifndef _BASICONEDMESH_H_
#define _BASICONEDMESH_H_


#include <cmath>
#include <vector>
#include "lifeV.hpp"
#include "GetPot.hpp"

namespace LifeV
{
/*!
  \class Point1D: geometrical point in 1D

*/
class Point1D
{
public:
  //! Default constructor
  Point1D();

  //! Constructor
  Point1D(const double& x, const UInt& id=0);

  //! Copy constructor
  Point1D(const Point1D& pt);

  //! operator=
  Point1D & operator= (const Point1D & pt);

  //! return the identity of the point
  INLINE UInt id() const {return _M_id;};

  //! return the abscissae (and other coordinates :=0)
  INLINE double x() const {return _M_x;};
  INLINE double y() const {return 0.;};
  INLINE double z() const {return 0.;};

protected:
  //! abscissa
  double _M_x;

  //! identity
  UInt _M_id;
};

/*!
  \class Edge1D: geometrical edge in 1D

*/
class Edge1D
{
public:
  //! Default constructor
  Edge1D::Edge1D();

  //! Constructor with two end abscissae
  Edge1D( const double& x1, const double& x2, const UInt& id=0 );

  //! Constructor with two end points
  Edge1D( const Point1D& pt1, const Point1D& pt2, const UInt& id=0 );

  //! operator=
  Edge1D & operator= (const Edge1D & edg);


  //! return the identity of the edge
  INLINE UInt id() const {return _M_id;};

  //! return the length of the edge
  INLINE double length() const {return _M_length;};

  //! return the end points
  Point1D pt1() const {return _M_pt1;};
  Point1D pt2() const {return _M_pt2;};

  // return one of the end points (i=1 or 2)
  Point1D point( int i ) const;

protected:
  //! first end point
  Point1D _M_pt1;
  //! second end point (pt1.x() < pt2.x())
  Point1D _M_pt2;

  //! edge length (always positive)
  double _M_length;

  //! identity
  UInt _M_id;
};

/*!
  \class BasicOneDMesh : very basic one d mesh handler.
         I don't have the courage to modify (and then debug!!)
	 the huge 3D Mesh class...
	 To be done...

*/
class BasicOneDMesh
{
public:
  //! constructor reading a mesh file: do it!
  BasicOneDMesh(const std::string& mesh_file,
		const std::string& mesh_dir);

  //! constructor for regular meshes
  BasicOneDMesh(const double& xl, const double& xr,
		const int& nx);

  //! return one edge of the list (BEWARE: start at 1)
  Edge1D edgeList( UInt iedg );

  //! return the full list of points
  const std::vector< Point1D >& pointList() const {return _M_pointList;};

  //! number of points
  UInt numVertices() const {return _M_pointList.size();};

  //! number of edges
  UInt numEdges() const {return _M_edgeList.size();};

  //! Output
  void showMe(std::ostream& c=std::cout, UInt verbose = 0);

protected:

  std::vector< Point1D >  _M_pointList;
  std::vector< Edge1D >  _M_edgeList;

};
}
#endif
