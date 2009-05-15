/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

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
#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>

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
    Point1D(const Real& x, const UInt& id=0);

    //! Copy constructor
    Point1D(const Point1D& pt);

    //! operator=
    Point1D & operator= (const Point1D & pt);

    //! return the identity of the point
    UInt id() const {return _M_id;}

    //! return the abscissae (and other coordinates :=0)
    Real x() const {return _M_x;}
    Real y() const {return 0.;}
    Real z() const {return 0.;}

  protected:
    //! abscissa
    Real _M_x;

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
    Edge1D();

    //! Copy constructor
    Edge1D(const Edge1D & edg);

    //! Constructor with two end abscissae
    Edge1D( const Real& x1, const Real& x2, const UInt& id=0 );

    //! Constructor with two end points
    Edge1D( const Point1D& pt1, const Point1D& pt2, const UInt& id=0 );

    //! operator=
    Edge1D & operator= (const Edge1D & edg);


    //! return the identity of the edge
    UInt id() const {return _M_id;};
    //! return the identity of the edge
    UInt localId() const {return 0;};

    //! return the length of the edge
    Real length() const {return _M_length;};

    //! return the end points
    Point1D pt1() const {return _M_pt1;};
    Point1D pt2() const {return _M_pt2;};

    // return one of the end points (i=1 or 2)
    Point1D point( const UInt& i ) const;

  protected:
    //! first end point
    Point1D _M_pt1;
    //! second end point (pt1.x() < pt2.x())
    Point1D _M_pt2;

    //! edge length (always positive)
    Real _M_length;

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

    /*!
      Constructor for regular/adaptive meshes. 
      Mesh refines around a point of abscissa alpha, where
      the model features geometrical / physical discontinuities

      \param xl left abscissa
      \param xr right abscissa
      \param nb_elem number of mesh elements
      \param alpha abscissa of the point of discontinuity
      \param delta length of the discontinuity zone
      \param n exponent for the n-th order polynomial law for the variation of 1D model parameters
      \param min_deltax minimum length of an interval
      \param adaptive boolean flag (0/1) to create regular/adaptive meshes
    */
    BasicOneDMesh(const Real& xl, const Real& xr,
		  const UInt& nb_elem,
		  const Real& alpha=0, const Real& delta=0,
		  const Real& n=0, 
		  const Real& min_deltax=1, const bool& adaptive=false);

    //! return one edge of the list (BEWARE: start at 1)
    Edge1D edgeList( const UInt& iedg ) const;

    //! return the one point of the list (BEWARE: start at 1)
    Point1D pointList( const UInt& iedg ) const ;

    //! return the full list of points
    const std::vector< Point1D >& pointList() const {return _M_pointList;};

    //! number of points
    UInt numVertices() const {return _M_pointList.size();};

    //! number of edges
    UInt numEdges() const {return _M_edgeList.size();};

    //! Output
    void showMe(std::ostream& c=std::cout, const UInt& verbose = 0);

  protected:

    std::vector< Point1D >  _M_pointList;
    std::vector< Edge1D >  _M_edgeList;

  };
}
#endif
