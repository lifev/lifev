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
#include "darcyHandler.hpp"
#include "chrono.hpp"

using namespace std;

DarcyHandler::DarcyHandler(const GetPot& data_file):
  DataDarcy(data_file),
  DataAztec(data_file),
  nbCoor(nDimensions),
  geoMap(geoBilinearHexa),
  qr(quadRuleHexa8pt),  
  geoMapBd(geoMap.boundaryMap()),
  qrBd(quadRuleQuad4pt),
  refBdFE(feQuadQ0),
  refVFE(feHexaRT0),
  refPFE(feHexaQ0),
  refTPFE(feHexaRT0Hyb),
  vfe(refVFE,geoMap,qr),
  pfe(refPFE,geoMap,qr),
  feBd(refBdFE,geoMapBd,qrBd),
  vdof(refVFE),pdof(refPFE),
  tpdof(refTPFE),
  numFacesPerVolume(mesh.volumeList(1).numLocalFaces) /* we assume that all
  element have the same number of faces, so we just look at the first one */
{
  // read mesh
  readINRIAMeshFile(mesh,mesh_dir+"/"+mesh_file,1);
  //  mesh.check(true,true);
  mesh.updateElementEdges();
  mesh.updateElementFaces(true); /*the "true" flag is to build the faceList
				   of all faces (and not only the boundary) */
  if(verbose>2) mesh.showMe();
  // build dof
  vdof.update(mesh);
  pdof.update(mesh);
  tpdof.update(mesh);
  
  dimTPdof = tpdof.numTotalDof();
  dimPdof = pdof.numTotalDof();
  dimVdof = vdof.numTotalDof();

  if(verbose>0){
    cout << "Number of TP dof : " << dimTPdof << endl;
    cout << "Number of  P dof : " << dimPdof << endl;
    cout << "Number of  V dof : " << dimVdof << endl;
  }
  // define the boundary conditions
  switch(test_case){
  case 1:
    /*
      Neumann condition (on p) at inlet (ref 1) and outlet (ref 3)
      example of mesh: cylhexa.mesh
    */
    bc_fct1.setFunction(g1);
    bc_fct2.setFunction(g3);
    nb_bc = 2;
    bc.setNumber(nb_bc);
    bc.addBC("Inlet",      1, Natural,   Scalar, bc_fct1); 
    bc.addBC("Outlet",  3, Natural, Scalar, bc_fct2); 
    break;

  case 2:
    /*
      Dirichlet condition (on p) at inlet (ref 1) and outlet (ref 3)
      example of mesh: cylhexa.mesh
    */      
    bc_fct1.setFunction(g1);
    bc_fct2.setFunction(g3);
    nb_bc = 2;
    bc.setNumber(nb_bc);
    bc.addBC("Inlet",      1, Essential,   Scalar, bc_fct1); 
    bc.addBC("Outlet",  3, Essential, Scalar, bc_fct2); 
    break;

  case 3:
    /*
      Robin condition (dp/dq + alpha p = 1, alpha=1) at inlet (ref 1) 
      and Dirichlet (p=-1) at outlet (ref 3)
      example of mesh: cylhexa.mesh
    */      
    bc_fct_rob.setFunctions_Mixte(g1, mixte_coeff); //! Robin coeff = 1.
    bc_fct2.setFunction(g3);
    nb_bc = 2;
    bc.setNumber(nb_bc);
    bc.addBC("Inlet",   1,     Mixte, Scalar, bc_fct_rob); 
    bc.addBC("Outlet",  3, Essential, Scalar, bc_fct2); 
    break;

  case 33:
    /* 
      Analytical solution defined in user_fct

      example of mesh: hexahexa10x10x10.mesh

    */
    bc_fct1.setFunction(zero); //!< low   X
    bc_fct2.setFunction(zero); //!< high  X
    bc_fct3.setFunction(zero); //!< low   Y
    bc_fct4.setFunction(zero); //!< high  Y
    bc_fct5.setFunction(zero); //!< low   Z
    bc_fct6.setFunction(zero); //!< high  Z

    nb_bc = 6;
    bc.setNumber(nb_bc);
    
    bc.addBC("Analytical, real BC",    1, Essential,   Scalar, bc_fct1); 
    bc.addBC("Analytical, real BC",    2, Essential,   Scalar, bc_fct2);
    bc.addBC("Analytical, real BC",    3, Essential,   Scalar, bc_fct3);
    bc.addBC("Analytical, real BC",    4, Essential,   Scalar, bc_fct4);
    bc.addBC("Analytical, real BC",    5, Essential,   Scalar, bc_fct5);
    bc.addBC("Analytical, real BC",    6, Essential,   Scalar, bc_fct6);
    
    break;  

  default:
    ERROR_MSG("Unknown test case");
  }
  // update the dof with the b.c.
  bc.bdUpdate(mesh, feBd, tpdof);
  //  tpdof.bdUpdate(mesh, feBd, bc);
  // check the mesh after b.c.
  /*
    BE CAREFUL: calling 
    mesh.check(true,true);
    after updateElementFaces(true) erase faceElement !!!!
  */
  //
  if(verbose>2) tpdof.showMe();
  if(verbose>2) bc.showMe(true);
}


