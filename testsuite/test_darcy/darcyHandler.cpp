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
  default:
    cerr << "Unknown test case\n";
    exit(1);
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


