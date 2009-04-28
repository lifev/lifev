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
#ifndef _GEOMAP_H
#define _GEOMAP_H

#include <life/lifecore/life.hpp>
#include <life/lifefem/refEle.hpp>
/*!
  \file geoMap.h
  \brief Structure for the geometrical mapping
*/

namespace LifeV
{
/*!
  \class GeoMap
  \brief Structure for the geometrical mapping
  \author J.-F. Gerbeau
  \date 04/2002

  This class contains the geometrical transformation that maps the reference
  element on the current element, and its values on integration points

  \par How to add a new geometrical mapping

  The way is very similar to a reference finite element see refFE.h
*/

class GeoMap:
            public RefEle
{
    const GeoMap* _boundaryMap;
public:
    //! Constructor of a geo map
    /*!
      Constructor of a geo map. The arguments are:

      _name : the name of the f.e.

      _shape : the geometry belongs to enum ReferenceShapes {NONE, POINT, LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA}; (see basisElSh.h)

       _nbDof : the total number of d.o.f.

       _nbCoor : number of local coordinates

       phi : the static array containing the basis functions (defined in refEle.h)

       dPhi : the static array containing the derivatives of the basis functions (defined in refEle.h)

       d2Phi : the static array containing the second derivatives of the basis functions (defined in refEle.h)

       refCoor : the static array containing the coordinates of the nodes on the reference element (defined in refEle.h)

       sqr : a set of quadrature rule (defined in quadRule.cc)

       bdMap : a pointer on the natural associated mapping for the boundary of the element
     */
    GeoMap( std::string _name, ReferenceShapes _shape, int _nbDof, int _nbCoor,
            const Fct* phi, const Fct* dPhi, const Fct* d2Phi,
            const Real* _refCoor, const SetOfQuadRule& sqr, const GeoMap* bdMap );
    ~GeoMap();
    friend std::ostream& operator << ( std:: ostream& f, const GeoMap& geomap );
    //! return the natural mapping for the boundary of the element
    inline const GeoMap& boundaryMap() const
    {
        ASSERT_PRE( _boundaryMap , "No boundary map defined" );
        return *_boundaryMap;
    }
};

//--------------------------------------------------
extern const GeoMap geoLinearSeg;
extern const GeoMap geoLinearTria;
extern const GeoMap geoBilinearQuad;
extern const GeoMap geoLinearTetra;
extern const GeoMap geoBilinearHexa;
//
/*! Helper function that returns the geomap associated to a mesh

\note To be completed!
*/
template <typename RegionMesh>
const GeoMap& getGeoMap( RegionMesh & /*mesh*/ )
{
    typedef typename RegionMesh::ElementShape ElementShape;
    switch ( ElementShape::Shape )
    {
    case HEXA:
        if ( ElementShape::numPoints == 8 )
            return geoBilinearHexa;
        else
            ERROR_MSG( "Geomap type not yet implemented" );
        break;
    case TETRA:
        if ( ElementShape::numPoints == 4 )
            return geoLinearTetra;
        else
            ERROR_MSG( "Geomap type not yet implemented" );
        break;
    case TRIANGLE:
            if ( ElementShape::numPoints == 3 )
                return geoLinearTria;
            else
                ERROR_MSG( "Geomap type not yet implemented" );
            break;
    case QUAD:
            if ( ElementShape::numPoints == 4 )
                return geoBilinearQuad;
            else
                ERROR_MSG( "Geomap type not yet implemented" );
            break;
    default:
        ERROR_MSG( "Geomap type not yet implemented" );
    }
}
}
#endif
