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
#ifndef _GEOMAPDG_H
#define _GEOMAPDG_H

#include <life/lifecore/life.hpp>
#include <life/lifefem/refEleDG.hpp>

/*!
  \file geoMapDG.h
  \brief Structure for the geometrical mapping
*/


namespace LifeV
{
/*!
  \class GeoMapDG
  \brief Structure for the geometrical mapping
  \author D. A. Di Pietro
  \date 12/2004

  This class contains the geometrical transformation that maps the reference
  discontinuos element on the current element, and its values on quadrature nodes.

  \par How to add a new geometrical mapping

  The way is very similar to a reference finite element see refFEDG.h
*/

class GeoMapDG:public RefEleDG{
 protected:
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
  GeoMapDG::GeoMapDG(std::string _name,
             ReferenceShapes _shape,
             int _nbDof,
             int _nbCoor,
             const Fct* phi, const Fct* dPhi, const Fct* d2Phi,
             const Real* refCoor,
             const SetOfQuadRule& sqr,
             ReferenceShapes _shapeFaces,
             int _nbFaces, int _nbGeoNodeFaces,
             const Real* refCoorFaces,
             const SetOfQuadRule& sqrFaces, const GeoMap& _geoMap, const GeoMap* bdMap);
  ~GeoMapDG();

  friend std::ostream& operator << (std::ostream& f,const GeoMapDG& geomap);

  //! return the natural mapping for the boundary of the element
  inline const GeoMap& boundaryMap() const
    {
      ASSERT_PRE( _boundaryMap , "No boundary map defined");
      return *_boundaryMap;
    }
};

//--------------------------------------------------
extern const GeoMap geoDGLinearSeg;
extern const GeoMap geoDGLinearTria;
extern const GeoMap geoDGBilinearQuad;
extern const GeoMap geoDGLinearTetra;
extern const GeoMap geoDGBilinearHexa;
#endif
}
