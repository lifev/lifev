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
  \file oneDModelHandler.cpp
  \author Vincent Martin
  \date 07/2004 
  \version 1.0

  \brief This file contains a very basic one d mesh handler.

*/

#include "basicOneDMesh.hpp"

//------------------------------------
//!  implementation for Point1D
//------------------------------------
//! Default constructor   
Point1D::Point1D():
  _M_x(0.),
  _M_id(0)
{ 
}
//! Constructor   
Point1D::Point1D(const double& x, const UInt& id):
  _M_x(x),
  _M_id(id)
{
}
//! Copy constructor   
Point1D::Point1D(const Point1D& pt):
  _M_x(pt._M_x),
  _M_id(pt._M_id)
{
}
//! operator =
Point1D & Point1D::operator= (const Point1D & pt)
{
  _M_x  = pt.x();
  _M_id = pt.id();
  return *this;
}

//------------------------------------
//!  implementation for Edge1D
//------------------------------------
//! Default constructor 
Edge1D::Edge1D():
  _M_pt1(),
  _M_pt2(),
  _M_length(0),
  _M_id(0)
{}

//! Constructor with two end abscissae
Edge1D::Edge1D( const double& x1, const double& x2, const UInt& id ):
  _M_pt1(std::min(x1,x2)),
  _M_pt2(std::max(x2,x1)),
  _M_length(std::abs(x2-x1)),
  _M_id(id)
{
  ASSERT_PRE(_M_length > 0 , 
	     "The points must be provided with (stricly) increasing absissa.");
}

//! Constructor with two end abscissae
Edge1D::Edge1D( const Point1D& pt1, const Point1D& pt2, const UInt& id ):
  _M_id(id)
{ 
  if ( double diff = pt2.x() - pt1.x() > 0 ) {
    _M_pt1 = pt1;  _M_pt2 = pt2;
    _M_length = diff;
  }
  else {
    _M_pt1 = pt2;  _M_pt2 = pt1;
    _M_length = -diff;   
  }
}

//! operator =
Edge1D & Edge1D::operator= (const Edge1D & edg)
{
  _M_pt1  = edg.pt1();
  _M_pt2  = edg.pt2();
  _M_length = edg.length();
  _M_id   = edg.id();
  return *this;
}

// return one of the end points (i=1 or 2)
Point1D Edge1D::point( int i ) const
{
  switch(i){
  case 1: 
    return _M_pt1; break;
  case 2: 
    return _M_pt2; break;
  default:
    ERROR_MSG("Only two end points for an edge (basic Edge1D)");
  }
}

//------------------------------------
//!  implementation for BasicOneDMesh
//------------------------------------
//! constructor reading a mesh file: do it!
BasicOneDMesh::BasicOneDMesh(const std::string& mesh_file,
			     const std::string& mesh_dir)
{}

//! constructor for regular meshes
BasicOneDMesh::BasicOneDMesh(const double& xl, const double& xr,
			     const int& nx):
  _M_pointList( nx + 2 ),
  _M_edgeList( nx + 1 )
{
  ASSERT_PRE( nx > 0, "The number of elements must be positive!");

  double x_current = xl;
  double deltax = ( xr - xl ) / (nx+1);
  ASSERT_PRE( deltax > 0 , 
	      "The left point is on the right..." );

  for (UInt it=0; it < _M_pointList.size(); it++)
    {
      _M_pointList[it] = Point1D(x_current, it);
      x_current += deltax;
    }
  ASSERT_PRE( std::abs(xr - x_current + deltax ) < 1e-10 * deltax , 
       "Some problems with the basic 1D mesh build. Check xleft<xright?" );

  for (UInt it=0; it < _M_edgeList.size(); it++)
    {
      _M_edgeList[it] = Edge1D(_M_pointList[it],_M_pointList[it+1], it);
    }
}

//! return one edge of the list (iedg starts at 1)
Edge1D BasicOneDMesh::edgeList( UInt iedg )
{ 
  ASSERT_BD(0 < iedg && iedg < _M_edgeList.size() + 1 );
  return _M_edgeList[ iedg - 1 ];
}

void BasicOneDMesh::showMe(std::ostream& c, UInt verbose)
{
  c << "\n*** Basic 1D Mesh \n";
  c << "number of points = " << numVertices() << "\n";
  if (verbose > 3) {
    UInt count(0),lines(10);
    for (UInt it=0; it < _M_pointList.size(); it++) {
      if (count++ % lines ==0){
	c << std::endl;
      }
      c << _M_pointList[it].x() << " " << _M_pointList[it].id() << "\t" ;
    }
  }
  c << "\nnumber of edges = " << numEdges() << "\n";
  if (verbose > 3) {
    for (UInt it=0; it < _M_edgeList.size(); it++) {
      c << "x1 = " << _M_edgeList[it].pt1().x()
	<< "\tx2 = " << _M_edgeList[it].pt2().x()
	<< "\t" << _M_edgeList[it].id()  << "\n";
    }
      
  }
  c << "*** End of Basic 1D Mesh \n" << std::endl;
}
