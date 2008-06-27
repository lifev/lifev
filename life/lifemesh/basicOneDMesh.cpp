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
  \file basicOneDMesh.cpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief This file contains a very basic one d mesh handler.

*/

#include <life/lifemesh/basicOneDMesh.hpp>

namespace LifeV
{
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
  Point1D::Point1D(const Real& x, const UInt& id):
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

  //! Copy constructor
  Edge1D::Edge1D(const Edge1D & edg) :
    _M_pt1( edg._M_pt1 ),
    _M_pt2( edg._M_pt2 ),
    _M_length( edg._M_length ),
    _M_id( edg._M_id )
  {}

  //! Constructor with two end abscissae
  Edge1D::Edge1D( const Real& x1, const Real& x2, const UInt& id ):
    _M_pt1(std::min(x1,x2)),
    _M_pt2(std::max(x2,x1)),
    _M_length(std::abs(x2-x1)),
    _M_id(id)
  {
    ASSERT_PRE(_M_length > 0 ,
               "BasicOneDMesh::BasicOneDMesh problems! The points must be provided with (stricly) increasing absissa.");
  }

  //! Constructor with two end abscissae
  Edge1D::Edge1D( const Point1D& pt1, const Point1D& pt2, const UInt& id ):
    _M_id(id)
  {
    if ( Real diff = pt2.x() - pt1.x() > 0 ) {
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
  Point1D Edge1D::point( const UInt& i ) const
  {
    switch(i){
    case 1:
      return _M_pt1; break;
    case 2:
      return _M_pt2; break;
    default:
      ERROR_MSG("Only two end points for an edge (basic Edge1D)");
    }
    return Point1D();
  }

  //------------------------------------
  //!  implementation for BasicOneDMesh
  //------------------------------------
  //! constructor reading a mesh file: do it!
  BasicOneDMesh::BasicOneDMesh(const std::string& /*mesh_file*/,
			       const std::string& /*mesh_dir*/)
  {}

  //! constructor for regular/adaptive mesh TP 1/05
  BasicOneDMesh::BasicOneDMesh(const Real& xl, const Real& xr,
			       const UInt& nb_elem,			     
			       const Real& alpha, const Real& delta,
			       const Real& n, 
			       const Real& min_deltax, const bool& adaptive):
    _M_pointList( nb_elem + 1 ),
    _M_edgeList( nb_elem )
  {
    if (adaptive) // algorithm for adaptive mesh
      {
	ASSERT_PRE( nb_elem > 0, "BasicOneDMesh::BasicOneDMesh problems! The number of elements must be positive!");
	ASSERT_PRE( ( ( alpha + delta/2 ) > xr ) || 
		    ( ( alpha - delta/2 ) < xl ) , "BasicOneDMesh::BasicOneDMesh problems! Invalid alpha and delta values!");
      
	UInt it=0, it_l, alpha_it;
	/*
	  we want to place as few points as possible outside the interval
	  ( alpha-delta/2, alpha+delta/2 )

	  |---------------|------------|------------|---------------|
	  xl         alpha-delta/2   alpha     alpha+delta/2        xr

	  given the size of point list (nb_elem + 1) we expect to have
	  -> [ (alpha-delta/2) / min_deltax ] points in (xl, alpha-delta/2)
	  -> [ (xr - (alpha+delta/2)) / min_deltax ] points in (alpha+delta/2, xr)
	  -> all the remaining points in (alpha-delta/2, alpha+delta/2)
	*/      
	// label of the point of abscissa alpha
	alpha_it = static_cast<int>( std::floor( 0.5 + (alpha-delta/2) / min_deltax ) ) +
	  ( _M_edgeList.size() - 
	    static_cast<int>( std::floor( 0.5 + (xr - (alpha+delta/2)) / min_deltax ) ) -
	    static_cast<int>( std::floor( 0.5 + (alpha-delta/2) / min_deltax ) ) ) / 2;
	// auxiliary variables
	Real ratio, n_elem_delta, n_elem_r,
	  x_current, x_current_l,
	  deltax, deltax_adaptive, deltax_uniform;
      	// number of mesh elements in the interval ( (alpha+delta/2) - alpha )
	n_elem_r = ( _M_edgeList.size() - alpha_it ) - 
	  static_cast<int>( std::floor( 0.5 + (xr - (alpha+delta/2)) / min_deltax ) );
	/*
	  number of mesh elements in the interval ( (alpha+delta/2) - (alpha-delta/2) )
	  when all elements have the same length
	*/
	n_elem_delta = static_cast<Real>( _M_edgeList.size() ) / (xr - xl) * delta;
      
	// Start from discontinuity point
	x_current = alpha;
	do
	  {
	    /*
	      ratio stands for the "percentage" of interval (alpha-delta/2, alpha+delta/2)
	      covered when considering x_current
	    */
	    ratio=(( (alpha + delta/2) - x_current ) / delta);
	    // insert a new Point1D in point list	  
	    _M_pointList[alpha_it+it] = Point1D(x_current, alpha_it+it);
	    // insert the symmetrical (with respect to alpha) Point1D in point list
	    _M_pointList[alpha_it-it] = Point1D( ( 2*alpha - x_current ),
						 alpha_it-it);
	    // evaluate interval size
	    /*
	      we expect to have one (or more) discontinuous physical / geometrical
	      property P: alpha is the discontinuity point.

	      We approximate the discontinuity with a smooth function f(x) \in (0, 1)
	      such that:
	      -> in (xl, alpha-delta\2), P = P_1
	      -> in (alpha+delta\2, xr), P = P_2
	      -> in (alpha-delta\2, alpha+delta/2), P = P_1 + f(x) * (P_2 - P_1)

	      Adaptive interval size (dx) is evaluated by imposing that
	      -> f(x + dx) - f(x) = 1 / n_elem_delta
	      By setting f(x + dx) = f(x) + f'(x) dx we find
	      -> dx = ( 1 / n_elem_delta )/ f'(x)

	      In this implementation we take
	      -> f(x) = 1 - 2^(n-1) * {[x - (alpha-delta/2)] * 1/delta }^n
	              for x in (alpha-delta\2, alpha)
	      -> f(x) = 2^(n-1) * {[(alpha+delta/2) - x] * 1/delta }^n
	              for x in (alpha, alpha+delta\2)

	      and, consequently:
	      -> f'(x) = - (n/delta) * 2^(n-1) * {[x - (alpha-delta/2)] * 1/delta}^(n-1)
	              for x in (alpha-delta\2, alpha)
	      -> f'(x) = - (n/delta) * 2^(n-1) * {[(alpha+delta/2) - x] * 1/delta}^(n-1)
	              for x in (alpha, alpha+delta\2)

	      Note that f(x) is symmetrical in (alpha-delta\2, alpha+delta/2), with respect
	      to point alpha. So we only compute deltax_adaptive inside (alpha-delta\2, alpha)
	    */
	    deltax_adaptive = ( - 1/n_elem_delta ) * 
	      ( 1 / ( (-n/delta) * 
		      pow( 2 * ratio , (n-1) ) 
		      ) 
		);
	    // uniform interval size
	    deltax_uniform = ( (alpha+delta/2) - x_current) / ( n_elem_r - it );
	    // point label
	    it++;
	    /*
	      this is made in order not to have too few intervals:
	      deltax_adaptive tends to high values when f'(x) and this could cause the mesh
	      to have less than n_elem_delta intervals inside (alpha-delta\2, alpha+delta/2).
	      So, when (deltax_adaptive > deltax_uniform) resort to deltax_uniform
	    */
	    deltax = ( (deltax_adaptive < deltax_uniform) && (it < n_elem_r) ) 
	      ? deltax_adaptive : deltax_uniform;
	  
	    ASSERT_PRE( deltax > 0 ,
			"BasicOneDMesh::BasicOneDMesh problems! The left point is on the right..." );
	    // move to the next point location
	    x_current += deltax;
	  }
	while ( x_current < ( alpha + delta/2 ) );
      
	it_l=it;
	x_current_l=x_current;
      
	if ( ( x_current ) < xr ) // we are in interval ( alpha + delta/2, xr )
	  {
	    // here the mesh is regular
	    do
	      {
		// insert a new Point1D in point list	  
		_M_pointList[alpha_it+it] = Point1D(x_current, alpha_it+it);
		/*
		  uniform interval size:
		  subdivide ( xr - x_current ) into ( _M_edgeList.size() - (alpha_it + it) )
		  intervals of the same size
		 */
		deltax = ( xr - x_current ) / ( _M_edgeList.size() - 
						(alpha_it + it) );
		// move to the next point location
		x_current += deltax;
		// update point label
		it++;
	      }
	    while ( ( x_current < xr ) && ( alpha_it+it < nb_elem ) );
	  
	    ASSERT_PRE( std::abs(xr - x_current + deltax ) < 1e-10 * deltax ,
			"BasicOneDMesh::BasicOneDMesh problems! Check xleft<xright?" );
	  }
	// add last point (boundary point xr)
	_M_pointList[nb_elem] = Point1D(xr, nb_elem);
      
	/*
	  at this level, x_current_l >= ( alpha + delta/2 ) so that 
	  -> 2*alpha - x_current_l <= alpha - delta/2
	 */
	if ( ( 2*alpha - x_current_l ) > xl ) // we are in interval ( xl, alpha - delta/2 )
	  {
	    // here the mesh is regular
	    do
	      {
		// insert a new Point1D in point list	  
		_M_pointList[alpha_it-it_l] = Point1D( ( 2*alpha - x_current_l ),
						       alpha_it-it_l);
		/*
		  uniform interval size:
		  subdivide ( (2*alpha - x_current_l ) - xl ) into ( alpha_it - it_l )
		  intervals of the same size
		 */
		deltax = ( ( 2*alpha - x_current_l ) - xl ) / ( alpha_it - it_l );
		// move to the next point location (towards xl)
		x_current_l += deltax;
		// update point label
		it_l++;
	      }
	    while ( (2*alpha - x_current_l) > 0 ) ;
	  
	    ASSERT_PRE( std::abs(xr - x_current + deltax ) < 1e-10 * deltax ,
			"BasicOneDMesh::BasicOneDMesh problems! Check xleft<xright?" );
	  }
	// add first point (boundary point xl)
	_M_pointList[0] = Point1D(xl, 0);
 	// store edge list    
	for (UInt it=0; it < _M_edgeList.size(); it++)
	  {
	    _M_edgeList[it] = Edge1D(_M_pointList[it],_M_pointList[it+1], it);
	  }
      }
    else // algorithm for regular mesh
      {
	ASSERT_PRE( nb_elem > 0, "BasicOneDMesh::BasicOneDMesh problems! The number of elements must be positive!");
	// start from boundary point xl
	Real x_current = xl;
	/*
	  uniform interval size:
	  subdivide ( xr - x_l ) into ( _M_edgeList.size() ) intervals of the same size
	*/
	Real deltax = ( xr - xl ) / _M_edgeList.size();
	ASSERT_PRE( deltax > 0 ,
		    "BasicOneDMesh::BasicOneDMesh problems! The left point is on the right..." );
      
	for (UInt it=0; it < _M_pointList.size(); it++)
	  {
	    // insert a new Point1D in point list	  
	    _M_pointList[it] = Point1D(x_current, it);
	    // move to the next point location
	    x_current += deltax;
	  }
	ASSERT_PRE( std::abs(xr - x_current + deltax ) < 1e-10 * deltax ,
		    "BasicOneDMesh::BasicOneDMesh problems! Check xleft<xright?" );
	
	// store edge list    
	for (UInt it=0; it < _M_edgeList.size(); it++)
	  {
	    _M_edgeList[it] = Edge1D(_M_pointList[it],_M_pointList[it+1], it);
	  }
      }
  }


  //! return one edge of the list (iedg starts at 1)
  Edge1D BasicOneDMesh::edgeList( const UInt& iedg ) const
  {
    ASSERT_BD(0 < iedg && iedg < _M_edgeList.size() + 1 );
    return _M_edgeList[ iedg - 1 ];
  }

  //! return the one point of the list (BEWARE: start at 1)
  Point1D BasicOneDMesh::pointList( const UInt& ipt ) const 
  {
    ASSERT_BD(0 < ipt && ipt < _M_pointList.size() + 1 );
    return _M_pointList[ ipt - 1 ];
  }

  void BasicOneDMesh::showMe(std::ostream& c, const UInt& verbose)
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
}
