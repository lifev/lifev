/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/*!
  \file oneDModelHandler.cpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief This file contains the implementation for the OneD model handler.

*/

#include "oneDModelHandler.hpp"

namespace LifeV
{
//! Constructor
OneDModelHandler::OneDModelHandler(const GetPot& data_file):
     DataOneDModel(data_file),
     _M_nbCoor(1),
     _M_mesh(_M_x_left,_M_x_right,_M_nb_elem),
     _M_geoMap(geoLinearSeg),
     _M_qr(quadRuleSeg3pt),
     _M_refFE(feSegP1),
     _M_dof1D(_M_mesh.numVertices()),
     _M_dimDof(_M_dof1D.numTotalDof()),
     _M_fe(_M_refFE,_M_geoMap,_M_qr),
     _M_GracePlot( data_file )
{
  /* Useless as long as we don't have a 1d mesh reader and handler...
  //! read mesh
  readINRIAMeshFile(_M_mesh,_M_mesh_dir+"/"+_M_mesh_file,1);
  //  _M_mesh.check(true,true);
  _M_mesh.updateElementEdges();

  _M_dimDof = _M_dof1D.numTotalDof();
  */

  std::cout << "dimDof = "  << _M_dimDof << std::endl;

}

void OneDModelHandler::showMeHandler(std::ostream& c, UInt verbose)
{
  c << "\n--- One Dimensional Model Handler\n";
  _M_mesh.showMe(c, verbose);
  c << "--- End of One Dimensional Model Handler\n";

}

}
