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

  \brief This file contains the implementation for the OneD model solver.

*/

#include "oneDModelHandler.hpp"

namespace LifeV
{
//! Constructor
OneDModelHandler::OneDModelHandler(const GetPot& data_file):
     DataOneDModel(data_file),
     DataAztec(data_file),
     _M_nbCoor(1),
     _M_mesh(_M_x_left,_M_x_right,_M_nx),
     _M_geoMap(geoLinearSeg),
     _M_qr(quadRuleSeg3pt),
     _M_refFE(feSegP1),
     _M_dof(_M_mesh.numVertices()),
     _M_dimDof(_M_mesh.numVertices()),
     _M_fe(_M_refFE,_M_geoMap,_M_qr),
     _M_GracePlot()
{
  /* Useless as long as we don't have a 1d mesh reader and handler...
  //! read mesh
  readINRIAMeshFile(_M_mesh,_M_mesh_dir+"/"+_M_mesh_file,1);
  //  _M_mesh.check(true,true);
  _M_mesh.updateElementEdges();

  _M_dimDof = _M_dof.numTotalDof();
  */

  std::cout << "dimDof = "  << _M_dimDof << std::endl;

}


void OneDModelHandler::showMeHandler(std::ostream& c, UInt verbose)
{
  c << "\n--- One Dimensional Model Handler\n";
  _M_mesh.showMe(c, verbose);
  c << "--- End of One Dimensional Model Handler\n";

}

//! Returns the concentration
ScalUnknown<Vector>& OneDModelHandler::AreaUnkn()
{
  return _M_AreaUnkn;
}

// ! Initialize with constant initial conditions concentration
void OneDModelHandler::initialize(const double& u0)
{
  _M_U_initial = u0;
}

// ! Initialize when  initial conditions concentration
void OneDModelHandler::initialize(const Function& c0, Real t0, Real dt)
{

}

// ! Initialize when initial values for the concentration are read from file
void OneDModelHandler::initialize(const std::string & vname)
{

  std::fstream Resfile(vname.c_str(),ios::in | ios::binary);
  if (Resfile.fail()) {std::cerr<<" Error in initialize: File not found or locked"<<std::endl; abort();}
  Resfile.read((char*)&_M_AreaUnkn(1),_M_AreaUnkn.size()*sizeof(double));
  Resfile.close();
}

}
