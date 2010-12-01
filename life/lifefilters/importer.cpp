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
   \file importMesh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-11-06
 */
#include <life/lifefilters/readMesh3D.hpp>
#include <life/lifefilters/readMesh2D.hpp>
#include <life/lifefilters/importer.hpp>


namespace LifeV
{

namespace detail
{
template<typename Elt>
void
import( std::string const& __filename, MeshFormat const& __format,
        RegionMesh3D<Elt> & mesh, EntityFlag regionFlag )
{
    switch ( __format )
    {
    case MESHPP:
        readMppFile( mesh, __filename, regionFlag );
        break;

    case INRIA:
        readINRIAMeshFile( mesh, __filename, regionFlag );
        break;

    case GMSH:
        readGmshFile( mesh, __filename, regionFlag );
        break;

    case NETGEN:
        readNetgenMesh( mesh, __filename, regionFlag );
        break;

    }

}
template<typename Elt>
void
import( std::string const& __filename, MeshFormat const& __format,
        RegionMesh2D<Elt> & mesh, EntityFlag regionFlag )
{
    switch ( __format )
    {
    case MESHPP:
    case INRIA:
    case NETGEN:
    {
        std::ostringstream __ostr;
        __ostr << "Unsupported file format for RegionMesh2D";
        throw std::invalid_argument( __ostr.str() );
    }
    break;
    case GMSH:
        readGmshFile( mesh, __filename, regionFlag );
        break;
    }

}
}
void
importer::import( RegionMesh2D<LinearTriangle> & mesh, EntityFlag regionFlag )
{
    detail::import( _M_filename, _M_format, mesh, regionFlag );
}
void
importer::import( RegionMesh2D<LinearQuad> & mesh, EntityFlag regionFlag )
{
    detail::import( _M_filename, _M_format, mesh, regionFlag );
}
void
importer::import( RegionMesh3D<LinearTetra> & mesh, EntityFlag regionFlag )
{
    detail::import( _M_filename, _M_format, mesh, regionFlag );
}
void
importer::import( RegionMesh3D<LinearHexa> & mesh, EntityFlag regionFlag )
{
    detail::import( _M_filename, _M_format, mesh, regionFlag );
}
}
