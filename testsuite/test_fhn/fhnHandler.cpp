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
#include <iostream>
#include "fhnHandler.hpp"
#include "chrono.hpp"

namespace LifeV
{

FhNHandler::FhNHandler(const GetPot& data_file):
  DataFhN(data_file),
  DataAztec(data_file),
  DataTransient(data_file),
  nbCoor(nDimensions),
  geoMap(geoLinearTetra),
  qr(quadRuleTetra15pt),
  refFE(feTetraP1),
  fe(refFE,geoMap,qr),
  geoMapBd(geoMap.boundaryMap()),
  qrBd(quadRuleTria3pt),
  refBdFE(feTriaP1),
  feBd(refBdFE,geoMapBd,qrBd),
  dof(refFE)
{
  // read mesh
  readINRIAMeshFile(mesh,mesh_dir+"/"+mesh_file,1);
  
  if(verbose>2) mesh.showMe();
  // build dof
  dof.update(mesh);

  dimdof = dof.numTotalDof();

  if(verbose>0){
    std::cout << "Number of dof : " << dimdof << std::endl;
  }
  // define the boundary conditions
  switch(test_case){
  case 1:
    {
    /*
      Dirichlet (or Robin) BC on reference 1 (inlet)
      example of mesh : recduct.mesh
    */
    nb_bc = 1;
    bc.setNumber(nb_bc);
    vector<ID> comp(1);
    comp[0] = 1;
    // coef u + diff du/dn = g 
    bc_fct_rob.setFunctions_Mixte(one,zero); // first argument: g , second : coef
    bc.addBC("In",1, Mixte,Component,bc_fct_rob,comp);
    break;
    }
  case 2:
    {
      /*
	example of mesh : coeur1959.mesh
	Homogeneous Neumann b.c.
      */
      break;
    }
  case 3:
    {
      /*
	example of mesh : coeur1959.mesh
	Robin b.c. to take into account the position of an electrode
      */
      nb_bc = 1;
      bc.setNumber(nb_bc);
      vector<ID> comp(1);
      comp[0] = 1;
      // coef u + diff du/dn = g 
      bc_fct_rob.setFunctions_Mixte(stim_g,stim_coef);
      // first argument: g , second : coef (see user_fct.cpp)
      bc.addBC("In",0, Mixte,Component,bc_fct_rob,comp);
      break;
    }
  default:
    ERROR_MSG("Unknown test case");
  }
  // update the dof with the b.c.
  bc.bdUpdate(mesh, feBd, dof);
  if(verbose>2) dof.showMe();
  if(verbose>2) bc.showMe(true);
}
}

