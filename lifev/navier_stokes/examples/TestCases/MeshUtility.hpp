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

#ifndef MESHUTILITY_HPP
#define MESHUTILITY_HPP

#include <string>
#include <iostream>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>


namespace LifeV
{

namespace MeshUtility
{

//! Read and partitioned a *.mesh file
void fillWithFullMesh ( boost::shared_ptr< RegionMesh<LinearTetra> >& meshPart,
                        const std::string& meshName,
                        const std::string& ressourcesPath = "./Ressources/" );

//! Read a partitioned mesh
void fillWithPartitionedMesh ( boost::shared_ptr< RegionMesh<LinearTetra> >& meshPart,
                               const std::string& meshName,
                               const std::string& ressourcesPath = "./Ressources/" );

//! Build a mesh from a partitioned mesh
void fillWithStructuredMesh ( boost::shared_ptr< RegionMesh<LinearTetra> >& meshPart,
                              markerID_Type regionFlag,
                              const UInt& m_x,
                              const UInt& m_y,
                              const UInt& m_z,
                              bool verbose = false,
                              const Real& l_x = 1.0,
                              const Real& l_y = 1.0,
                              const Real& l_z = 1.0,
                              const Real& t_x = 0.0,
                              const Real& t_y = 0.0,
                              const Real& t_z = 0.0 );

//! Print informations about the mesh
void printMeshInfos ( boost::shared_ptr<RegionMesh<LinearTetra> > mesh );

} // namespace MeshUtility

} // namespace LifeV

#endif /* PDEPROBLEM_HPP */
