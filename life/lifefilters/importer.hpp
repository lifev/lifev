/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-11-06

  Copyright (C) 2004 EPFL

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
/**
   \file importMesh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-11-06
 */
#ifndef __importer_H
#define __importer_H 1

#include <life/lifemesh/regionMesh3D.hpp>


namespace LifeV
{
enum MeshFormat {
    MESHPP,
    INRIA,
    GMSH,
    NETGEN
};
class importer
{
public:

    importer()
        :
        _M_filename(),
        _M_format(GMSH)
        {}
    importer( std::string const& filename,  MeshFormat const& format )
        :
        _M_filename( filename ),
        _M_format( format )
        {}

    void setFilename( std::string const& __filename ) {
        _M_filename = __filename;
    }

    void setFormat( MeshFormat const& __format ) {
        _M_format = __format;
    }

    //! import mesh with tetrahedras
    void import( RegionMesh3D<LinearTetra> & mesh, EntityFlag regionFlag );

    //! import mesh with hexahedras
    void import( RegionMesh3D<LinearHexa> & mesh, EntityFlag regionFlag );

private:

    //! name of the file to import
    std::string _M_filename;

    //! format of the file to import
    MeshFormat _M_format;
};
}

#endif /* __importer_H */
