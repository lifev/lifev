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
    @brief This file contains the definition of the GeoMap class (and an helper function)

    @author Jean-Frederic Gerbeau
    @date 00-04-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    This class contains the geometrical transformation that maps the reference
    element on the current element.
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


  Modified by S. Quinodoz (samuel.quinodoz@epfl.ch, 04.2010)
*/
class GeoMap:
        public RefEle
{
public:

    //! @name Public Types
    //@{

    typedef RefEle::function_Type function_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Full Constructor of a geo map
    /*!
      @param name : the name of the f.e.
      @param shape : the geometry belongs to enum ReferenceShapes {NONE, POINT, LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA}; (see basisElSh.h)
      @param nbDof : the total number of d.o.f.
      @param nbCoor : number of local coordinates
      @param phi : the static array containing the basis functions (defined in refEle.h)
      @param dPhi : the static array containing the derivatives of the basis functions (defined in refEle.h)
      @param d2Phi : the static array containing the second derivatives of the basis functions (defined in refEle.h)
      @param refCoor : the static array containing the coordinates of the nodes on the reference element (defined in refEle.h)
      @param  bdMap : a pointer on the natural associated mapping for the boundary of the element
     */
    GeoMap( std::string          name,
            ReferenceShapes      shape,
            UInt                 nbDof,
            UInt                 nbCoor,
            const function_Type* phi,
            const function_Type* dPhi,
            const function_Type* d2Phi,
            const Real*          refCoor,
            const GeoMap*        bdMap );

    //! Destructor
    ~GeoMap();

    //@}


    //! @name Get Methods
    //@{

    //! return the natural mapping for the boundary of the element
    const GeoMap& boundaryMap() const
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
template <typename MeshType>
const GeoMap& getGeoMap( MeshType & /*mesh*/ )
{

    typedef typename MeshType::ElementShape ElementShape;

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
