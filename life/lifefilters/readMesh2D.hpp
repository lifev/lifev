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
    @brief Mesh reader from mesh2d files

    @author Luca Formaggia <luca.formaggia@polimi.it>
    @contributor Nur Aiman Fadel <nur.fadel@mail.polimi.it>
    @maintainer Nur Aiman Fadel <nur.fadel@mail.polimi.it>

    @date 19-08-1999

    Mesh reader that it is able to read 2d meshes from freefem and gmsh files.<br>
    readMesh2d reads only triangle meshes.<br>
    readMppFiles handles only liner&quad triangles.<br>
 */

#ifndef _READMESH2D_HH_
#define _READMESH2D_HH_ 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/lambda/lambda.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/fortran_wrap.hpp>
#include <life/lifecore/util_string.hpp>

#include <life/lifemesh/regionMesh2D.hpp>

#include <life/lifefilters/mesh_util.hpp>

namespace LifeV
{

// ===================================================
// Macros for Fortran interface
// ===================================================

// Subroutine readmesh2d
SUBROUTINE_F77 F77NAME( readmesh2d ) ( I_F77 & ne, I_F77 & np,
                                       I_F77 & nptot, I_F77 & nb, I_F77 & nps, I_F77 & nx,
                                       I_F77 & ndimn, I_F77 & npe, I_F77 & npb,
                                       I_F77 & npc, I_F77 * iel,
                                       I_F77 & nd, R4_F77 * coor, I_F77 & ndc,
                                       R4_F77 & xmin, R4_F77 & xmax, R4_F77 & ymin,
                                       R4_F77 & ymax, I_F77 * ib, I_F77 & nbd,
                                       I_F77 * ic, I_F77 * bc, I_F77 * ie,
                                       I_F77 * cpl, R4_F77 * xmed, I_F77 & isw,
                                       I_F77 & ierr, CHARACTER filename );

// Subroutine read_mesh2d_head(filename,ne,np,nptot,npe,nb,nx,npc,ierr)
SUBROUTINE_F77 F77NAME( readmesh2dhead ) ( I_F77 & ne, I_F77 & np,
                                           I_F77 & nptot, I_F77 & npe, I_F77 & nb,
                                           I_F77 & nps, I_F77 & nx,
                                           I_F77 & npc, I_F77 & ierr, CHARACTER filename );

//! readMesh2d - reads a mesh in mesh2d (LF) format.
/*!
  It reads a gmsh mesh (2D) file and store it in a RegionMesh2D.

  @param mesh, the mesh data structure to fill in.
  @param fname, the name of the mesh file  to read.
  @param regionFlag, the identifier for the region.
  @return true if everything went fine, false otherwise.
*/

template <typename RegionMesh2D>
bool
readMesh2d( RegionMesh2D & mesh, const std::string & fname, EntityFlag regionFlag );

//! readGmshFile - reads a mesh in GMSH 2d format.
/*!
  It reads a gmsh mesh (2D) file and store it in a RegionMesh2D.

  @param mesh, the mesh data structure to fill in.
  @param filename, the name of the gmsh mesh file  to read.
  @param regionFlag, the identifier for the region.
  @return true if everything went fine, false otherwise.
*/

template <typename GeoShape, typename MC>

bool
readGmshFile( RegionMesh2D<GeoShape, MC> & mesh,
              const std::string & filename,
              EntityFlag regionFlag );

//! readFreeFemFile - reads a mesh in FreeFem 2d format.
/*!
read a freefem mesh (2D) file and store it in a RegionMesh2D.

@param mesh, the mesh data structure to fill in.
@param filename, the name of the freefem mesh file to read.
@param regionFlag, the identifier for the region.
@param bool useless, it will be removed.
@return true if everything went fine, false otherwise.
 */

template <typename GeoShape, typename MC>
bool
readFreeFemFile( RegionMesh2D<GeoShape, MC> & mesh,
                 const std::string & filename,
                 EntityFlag regionFlag,
                 bool useless );

} // Namespace LifeV

#endif /* READMESH2D_H */
