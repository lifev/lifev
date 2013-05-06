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
    @file MeshUtility.hpp
    @brief Functions to load meshes

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 16-01-2013
 */

#ifndef NSMESHUTILITY_HPP
#define NSMESHUTILITY_HPP

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <string>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/navier_stokes/examples/TestCases/MeshUtility.hpp>


namespace LifeV
{

namespace MeshUtility
{

//! setup and get the mesh data
/*!
  @param meshName name of the mesh file
  @param resourcesPath path to the mesh folder
  @param meshOrder order of the mesh elements
*/
/*!
    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
 */
MeshData
getMeshData ( const std::string& meshName,
              const std::string& resourcesPath,
              const std::string& meshOrder )
{
    MeshData meshData;
    meshData.setMeshDir ( resourcesPath );
    meshData.setMeshFile ( meshName );
    meshData.setMeshType ( ".mesh" );
    meshData.setMOrder ( meshOrder );
    meshData.setVerbose ( false );
    return meshData;
}

} // namespace MeshUtility

} // namespace LifeV

#endif /* NSMESHUTILITY_HPP */
