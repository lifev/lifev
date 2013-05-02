//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains an helper function to get a geometric map adapted to a mesh

     @date 07/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef MESH_GEOMETRIC_MAP_HPP
#define MESH_GEOMETRIC_MAP_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/GeometricMap.hpp>

#include <lifev/core/mesh/ElementShapes.hpp>

#include <boost/shared_ptr.hpp>


// LifeV namespace.
namespace LifeV
{


//! Generic implementation of the GeometricMapFromElementShape.
/*
This generic implementation should never be implemented. In case one asks
for a mapping that does not exist, an error will appear at compile time.

The specializations can be used to get the geometric mapping that goes
with the Shape given as template argument.
*/
template <typename Shape>
const GeometricMap& geometricMapFromElementShape();

//! \cond

template<>
inline const GeometricMap& geometricMapFromElementShape<GeoPoint>()
{
    return geoLinearNode;
}

template<>
inline const GeometricMap& geometricMapFromElementShape<LinearLine>()
{
    return geoLinearSeg;
}

template<>
inline const GeometricMap& geometricMapFromElementShape<LinearTriangle>()
{
    return geoLinearTria;
}

template<>
inline const GeometricMap& geometricMapFromElementShape<LinearQuad>()
{
    return geoBilinearQuad;
}

template<>
inline const GeometricMap& geometricMapFromElementShape<LinearTetra>()
{
    return geoLinearTetra;
}

template<>
inline const GeometricMap& geometricMapFromElementShape<LinearHexa>()
{
    return geoBilinearHexa;
}


//! \endcond

//! Function to get the map that goes with a mesh (version with template argument only)
/*!
    Given a type of mesh, this method returns an instance of geometric mapping
    that corresponds with the mesh.
*/
template <typename MeshType>
inline const GeometricMap& geometricMapFromMesh()
{
    return geometricMapFromElementShape<typename MeshType::elementShape_Type>();
}


//! Function to get the map that goes with a mesh (version with mesh in argument)
/*!
    Given a type of mesh, this method returns an instance of geometric mapping
    that corresponds with the mesh.
*/
template <typename MeshType>
inline const GeometricMap& geometricMapFromMesh (const boost::shared_ptr<MeshType>& /*mesh*/ )
{
    return geometricMapFromElementShape<typename MeshType::elementShape_Type>();
}

} // namespace LifeV


#endif

