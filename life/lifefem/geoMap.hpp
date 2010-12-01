//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief This file contains the definition of the GeoMap class (and an helper function)
 */

#ifndef GEOMAP_H
#define GEOMAP_H 1

#include <life/lifecore/life.hpp>

#include <life/lifefem/refEle.hpp>

namespace LifeV
{

//! GeoMap - Structure for the geometrical mapping
/*!
  @author J.-F. Gerbeau
  @date 04/2002

  This class contains the geometrical transformation that maps the reference
  element on the current element.

  Modified by S. Quinodoz (samuel.quinodoz@epfl.ch, 04.2010)
*/
class GeoMap:
        public RefEle
{
public:

    typedef RefEle::Fct Fct;

    //! @name Constructor & Destructor
    //@{

    //! Full Constructor of a geo map
    /*!
      @param _name : the name of the f.e.
      @param _shape : the geometry belongs to enum ReferenceShapes {NONE, POINT, LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA}; (see basisElSh.h)
      @param _nbDof : the total number of d.o.f.
      @param _nbCoor : number of local coordinates
      @param phi : the static array containing the basis functions (defined in refEle.h)
      @param dPhi : the static array containing the derivatives of the basis functions (defined in refEle.h)
      @param d2Phi : the static array containing the second derivatives of the basis functions (defined in refEle.h)
      @param refCoor : the static array containing the coordinates of the nodes on the reference element (defined in refEle.h)
      @param  bdMap : a pointer on the natural associated mapping for the boundary of the element
     */
    GeoMap( std::string          _name,
            ReferenceShapes      _shape,
            UInt                  _nbDof,
            UInt                  _nbCoor,
            const Fct*           phi,
            const Fct*           dPhi,
            const Fct*           d2Phi,
            const Real*          _refCoor,
            const GeoMap*        bdMap );

    //! Destructor
    ~GeoMap();

    //@}


    //! @name Get Methods
    //@{

    //! return the natural mapping for the boundary of the element
    inline const GeoMap& boundaryMap() const
    {
        ASSERT( M_boundaryMap!=0 , "No boundary map defined" );
        return *M_boundaryMap;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    GeoMap();

    //! No copy constructor
    GeoMap(const GeoMap&);

    //@}

    const GeoMap* M_boundaryMap;
};



//---- Predeclaration of the map (defined in defQuadRuleFE.cpp) ----

extern const GeoMap geoLinearNode;
extern const GeoMap geoLinearSeg;
extern const GeoMap geoLinearTria;
extern const GeoMap geoBilinearQuad;
extern const GeoMap geoLinearTetra;
extern const GeoMap geoBilinearHexa;


// ----

/*! Helper function that returns the geomap associated to a mesh */
template <typename RegionMesh>
const GeoMap& getGeoMap( RegionMesh & /*mesh*/ )
{

    typedef typename RegionMesh::ElementShape ElementShape;

    switch ( ElementShape::Shape )
    {
    case POINT:
        if ( ElementShape::numPoints == 1 )
            return geoLinearNode;
        else
            ERROR_MSG( "Geomap type not yet implemented" );
        break;
    case LINE:
        if ( ElementShape::numPoints == 2 )
            return geoLinearSeg;
        else
            ERROR_MSG( "Geomap type not yet implemented" );
        break;
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
