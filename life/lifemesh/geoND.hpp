/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*! file geoND.h */
#ifndef _GEOND_HH_
#define _GEOND_HH_

#include "geo0D.hpp"
#include "meshEntity.hpp"

namespace LifeV
{
/*!  Base class for Multi-dimensional basis Geometrical Entities.

  It has no boundary information, since the GeoXD boundary items will be
  ALWAYS stored first on the corresponding RegionMesh List. This is a
  paradigm of our code.

*/
template <typename GEOSHAPE,typename POINTTYPE=Geo0D>
class GeoND
    :
     public MeshEntity,
     public GEOSHAPE
{
public:

  GeoND();
  explicit GeoND(ID id);
  GeoND(const GeoND<GEOSHAPE,POINTTYPE> &);
  GeoND & operator=(GeoND const & G);
  ~GeoND();

  typedef GEOSHAPE GeoShape;
  typedef POINTTYPE  PointType;

  //! Number of points associated to the entity
  static const UInt numLocalPoints=GEOSHAPE::numPoints;// for comp only
  //! Number of Vertices associated to the entity
  static const UInt numLocalVertices=GEOSHAPE::numVertices;// for comp only
  //! The ith point (starting from 1)
  /* It returns the reference to an point object (possibly derived from
     Geo0D)*/
  PointType & point(ID const i);
  PointType const & point (ID const i) const;
  //! The ith point (starting from the end)
  /* It returns the reference to an point object (possibly derived from
  Geo0D). It starts from the last point, yet it follows the rule: vertices
  first. It may be used to access the points of a Geometry Element in a
  reverse way (i.e. with the opposite GeoElement orientation)*/
  PointType & reversepoint(ID const i);
  PointType const & reversepoint (ID const i) const;

  //!Inserts a point.  Uses point references
  void setPoint(ID const i,Geo0D  const & p); //put point
  //!Inserts a point Uses point references (bounds check)
  bool setPointBD(ID const i,Geo0D  const & p); //with forced bound check
  //!Inserts a point. Uses pointers
  void setPoint(ID const i,Geo0D  const * p); //put point
  //!Inserts a point. Uses pointers (bounds check)
  bool setPointBD(ID const i,Geo0D  const * p); //with forced bound check

  std::ostream &  showMe(bool verbose=false, std::ostream & c=std::cout) const;

  /*! Swap Points

    This is a member function to be used ONLY by routines for checking or
    amending meshes. You must give the local ID (which start from 1 as
    usual!)  */
  void swapPoints(const ID & pt1, const ID & pt2);

  /*! Exchange Points

    Exchanges points according to a list of old2new local ID numbering !
    old2new[i] is the new local ID of a point whose old local ID was ! i+1
    (remeber the numbering from 1 of the ID's!. This is a member function
    to be used ONLY by routines for checking or amending meshes. You must
    give ID (which start from 1 as usual!)  */
  void exchangePoints(const ID otn[GEOSHAPE::numPoints]);

 private:
  Geo0D *  _points[GEOSHAPE::numPoints];
};


/*--------------------------------------------------------------
                 GeoND
---------------------------------------------------------------*/
template <typename GEOSHAPE,typename POINTTYPE>
const UInt GeoND<GEOSHAPE,POINTTYPE>::numLocalPoints;

template <typename GEOSHAPE,typename POINTTYPE>
const UInt GeoND<GEOSHAPE,POINTTYPE>::numLocalVertices;


template <typename GEOSHAPE,typename POINTTYPE>
GeoND<GEOSHAPE,POINTTYPE>::GeoND():
MeshEntity(0)
{}

template <typename GEOSHAPE,typename POINTTYPE>
GeoND<GEOSHAPE,POINTTYPE>::GeoND(ID id):
MeshEntity(id)
{}

template <typename GEOSHAPE,typename POINTTYPE>
GeoND<GEOSHAPE,POINTTYPE>::GeoND(GeoND<GEOSHAPE,POINTTYPE> const & G):
  // Copy constructor
MeshEntity(G._id)
{
  for (UInt i =0; i< GeoND<GEOSHAPE,POINTTYPE>::numLocalPoints; ++i)
    _points[i]=G._points[i];
}

template <typename GEOSHAPE,typename POINTTYPE>
GeoND<GEOSHAPE,POINTTYPE>::~GeoND()
{}



template <typename GEOSHAPE,typename POINTTYPE>
GeoND<GEOSHAPE,POINTTYPE> &
GeoND<GEOSHAPE,POINTTYPE>::operator=(GeoND<GEOSHAPE,POINTTYPE> const & G)
     // Assignement operator
{
  if (this != &G){
    _id=G._id;
    for (UInt i =0; i< GeoND<GEOSHAPE,POINTTYPE>::numLocalPoints; ++i)
      _points[i]=G._points[i];
  }
  return *this;
}


template <typename GEOSHAPE,typename POINTTYPE>
INLINE
POINTTYPE & GeoND<GEOSHAPE,POINTTYPE>::point(ID const i)
{
  ASSERT_BD((i>0 && i<= GeoND<GEOSHAPE,POINTTYPE>::numLocalPoints)) ;
  return *(static_cast<POINTTYPE *>(_points[i-1])); // indexing from 1
};

template <typename GEOSHAPE,typename POINTTYPE>
INLINE
POINTTYPE const & GeoND<GEOSHAPE,POINTTYPE>::point(ID const i) const
{
  ASSERT_BD((i>0 && i<= GeoND<GEOSHAPE,POINTTYPE>::numLocalPoints));
  return *(static_cast<POINTTYPE *>(_points[i-1]));
};

template <typename GEOSHAPE,typename POINTTYPE>
INLINE
POINTTYPE & GeoND<GEOSHAPE,POINTTYPE>::reversepoint(ID const i)
{
  ASSERT_BD((i>0 && i<= GeoND<GEOSHAPE,POINTTYPE>::numLocalPoints)) ;
  return *(static_cast<POINTTYPE *>(_points[reversePoint<GEOSHAPE>::operate(i) -1])); // indexing from 1
};

template <typename GEOSHAPE,typename POINTTYPE>
INLINE
POINTTYPE const & GeoND<GEOSHAPE,POINTTYPE>::reversepoint(ID const i) const
{
  ASSERT_BD((i>0 && i<= GeoND<GEOSHAPE,POINTTYPE>::numLocalPoints));
  return *(static_cast<POINTTYPE *>(_points[reversePoint<GEOSHAPE>::operate(i) -1]));
};


template <typename GEOSHAPE,typename POINTTYPE>
INLINE
void  GeoND<GEOSHAPE,POINTTYPE>::setPoint(UInt const i,Geo0D  const & p)
{
  ASSERT_BD( (i > 0 && i <= GeoND<GEOSHAPE,POINTTYPE>::numLocalPoints )) ;
  _points[i-1]= const_cast<Geo0D  *>(&p);
};


template <typename GEOSHAPE,typename POINTTYPE>
bool  GeoND<GEOSHAPE,POINTTYPE>::setPointBD(UInt const i,Geo0D  const & p)
{
  ASSERT_BD0( (i > 0 && i <=  GeoND<GEOSHAPE,POINTTYPE>::numLocalPoints) ) ;

  // if not assert we need anyway to avoid under/overflows

  if ( i <= 0 || i > GeoND<GEOSHAPE,POINTTYPE>::numLocalVertices )
    return false;

  _points[i-1]=const_cast<Geo0D *>(&p);
  return true;
};

template <typename GEOSHAPE,typename POINTTYPE>
INLINE
void  GeoND<GEOSHAPE,POINTTYPE>::setPoint(UInt const i,Geo0D  const * p)
{
  ASSERT_BD( (i > 0 && i <= GeoND<GEOSHAPE,POINTTYPE>::numLocalPoints) ) ;
  _points[i-1]= const_cast<Geo0D  *>(p);
};


template <typename GEOSHAPE,typename POINTTYPE>
bool  GeoND<GEOSHAPE,POINTTYPE>::setPointBD(UInt const i,Geo0D  const * p)
{
  ASSERT_BD0( (i > 0 && i <=  GeoND<GEOSHAPE,POINTTYPE>::numLocalPoints) ) ;

  // if not assert we need anyway to avoid under/overflows

  if ( i <= 0 || i > GeoND<GEOSHAPE,POINTTYPE>::numLocalVertices )
    return false;

  _points[i-1]=const_cast<Geo0D *>(p);
  return true;
};


template <typename GEOSHAPE,typename POINTTYPE>
std::ostream &  GeoND<GEOSHAPE,POINTTYPE>::
showMe(bool verbose,std::ostream & out) const
{
  out << " GeoND object " <<std::endl;
  out << " Number of Vertices=" << GEOSHAPE::numVertices<<std::endl;
  out << " Number of Points=" << GEOSHAPE::numPoints<<std::endl;
  out << "ID= "<< id()<<std::endl;
  if (verbose)
    {
      out << " POINTS INFORMATION"<<std::endl<<std::endl;
      for (unsigned i=1 ; i<= GEOSHAPE::numVertices; i++)
	{
	  out << "POINT ID."<< i << std::endl;
	  out << point(i).showMe(verbose,out);
	}
    }
  out << "----- END OF GeoND data ---"<<std::endl<<std::endl;
  return out;
};

template <typename GEOSHAPE,typename POINTTYPE>
void  GeoND<GEOSHAPE,POINTTYPE>::swapPoints(const ID & pt1, const ID & pt2)
{
  Geo0D * tmp(_points[pt1-1]);
  _points[pt1-1]=_points[pt2-1];
  _points[pt2-1]=tmp;
}

template <typename GEOSHAPE,typename POINTTYPE>
void  GeoND<GEOSHAPE,POINTTYPE>::exchangePoints(const ID otn[GEOSHAPE::numPoints])
{
  Geo0D *  tmp[GEOSHAPE::numPoints];
  for (unsigned int i=0;i<GEOSHAPE::numPoints;++i){
    tmp[i]=_points[i];
  }
  for (unsigned int i=0;i<GEOSHAPE::numPoints;++i){
    _points[i]=tmp[otn[i]-1];
  }
}
}
#endif
