/*! --------------------------------------------------------------------------*
/                                                                            /
/                                                                            /
/ VTK WRITERS UTILITIES                                                      /
/                                                                            /
/ #Version 0.1 Experimental:                                                 /
/ 21/02/2001 Alessandro Veneziani                                            /
/ 10/07/2002 JFG                                                             /
/                                                                            /
/ #Purpose: Container for subroutines for VTK                                /
/                                                                            /
/                                                                            /
/---------------------------------------------------------------------------*/
#ifndef _VTK_WRTRS_
#define _VTK_WRTRS_

#define VTK_VERSION 2.3

// Linear cells
#define VTK_VERTEX 1
#define VTK_POLY_VERTEX 2
#define VTK_LINE 3
#define VTK_POLY_LINE 4
#define VTK_TRIANGLE 5
#define VTK_TRIANGLE_STRIP 6
#define VTK_POLYGON 7
#define VTK_PIXEL 8
#define VTK_QUAD 9
#define VTK_TETRA 10
#define VTK_VOXEL 11
#define VTK_HEXAHEDRON 12
#define VTK_WEDGE 13
#define VTK_PYRAMID 14

// Quadratic, isoparametric cells
#define VTK_QUADRATIC_EDGE       21
#define VTK_QUADRATIC_TRIANGLE   22
#define VTK_QUADRATIC_QUAD       23
#define VTK_QUADRATIC_TETRA      24
#define VTK_QUADRATIC_HEXAHEDRON 25

#include <fstream>
#include <vector>
#include "lifeV.hpp"
#include "currentFE.hpp"

/*!
  This subroutines writes the header part of a VTK file (which is constant in a time-dependent but fixed-mesh problem)
  for a general problem, in particular for instance for a P2 on LinearTetra Mesh: in such a case
  we need to compute the coordinates of all the supplementary nodes
  For "iso" cases (whenever we have ALL the coordinates of the nodes, it is recommended using the other subroutine
  with the same name but with a different signature.


  NB: it works with P1, P2, but it needs probably to be adapted to more general cases...
*/
template<typename TheMesh, typename TheDof>
void wr_vtk_ascii_header(string fname, string title,const  TheMesh& mesh, const TheDof& dof,CurrentFE& fe)
{
  // Note: it is assumed that all the (volume) element of the mesh are associated to the same reference finite element
 ofstream ofile(fname.c_str());

 const RefFE& refFE = fe.refFE;
 //  const GeoMap& geoMap = fe.geoMap;
 //  const QuadRule& qr = fe.qr;
 UInt cell_type;

 switch(refFE.type){
 case FE_P1_3D:
   cell_type = VTK_TETRA;
   break;
 case FE_P2_3D:case FE_P2tilde_3D:
   cell_type = VTK_QUADRATIC_TETRA; // maybe to be modified (P2 on linear tetra...) see also the geomap...
   break;
 case FE_Q1_3D:case FE_Q0_3D:
   cell_type = VTK_HEXAHEDRON; 
   break;
 default:
   cout << "WARNING: the element is not yet implemented in vtk_wrtrs.h\n";
   return;
 }
 
 ASSERT(ofile,"Error: Output file cannot be open"); //
 
 UInt nldpe=refFE.nbDofPerEdge;
 UInt nldpv=refFE.nbDofPerVertex;
 UInt nldpf=refFE.nbDofPerFace;
 UInt nldpV=refFE.nbDofPerVolume;
 UInt nlv=dof.numLocalVertices();
 UInt nle=dof.numLocalEdges();
 UInt nlf=dof.numLocalFaces();
 UInt nV=mesh.numVolumes();
 UInt ne=mesh.numEdges();
 UInt nv=mesh.numVertices();
 UInt nf=mesh.numFaces();
 UInt nldof = nle*nldpe+nlv*nldpv+nlf*nldpf+nldpV;
 //  _totalDof=nV*nldpV+ne*nldpe+nv*nldpv+nf*nldpf;
 UInt cells_size; 
 UInt num_points =  dof.numTotalDof(); // discuterne con Luca
 UInt num_points_supp =  dof.numTotalDof() - nv; // discuterne con Luca
 vector<Real> supp_x(num_points_supp,0.0);
 vector<Real> supp_y(num_points_supp,0.0);
 vector<Real> supp_z(num_points_supp,0.0);
 Real x,y,z;

 ofile << "# vtk DataFile Version " << VTK_VERSION << endl;
 ofile << title << endl;
 ofile << "ASCII" << endl;
 ofile << endl;
 ofile << "DATASET UNSTRUCTURED_GRID" << endl;
 ofile << "POINTS " << num_points  << " " << "float" << endl; // forse si puo' fare una RTTI sul dato contenuto nella point list
   
  UInt i,ie,index,j;
  UInt gcount;
  UInt lcount;

  // Vertex based Dof: the coordinates are available from the Pont List
  cout << "nv = " << nv << endl;
  for(i=0;i<nv;++i) // BUG ???????? i=1 ... !!!!!! (jfg 20/10/2002)
    ofile << mesh.pointList[i].x() << " "
	  <<  mesh.pointList[i].y() << " "
	  <<  mesh.pointList[i].z() << endl; 
 
 // Now I store the coordinates of the supplementary nodes in a temporary vector
  // Edge Based Dof
  gcount = 0;
  lcount = nlv;
  if (nldpe >0 ){
    for (ie=0; ie< nV; ++ie){
      fe.updateJac(mesh.volumeList(ie+1)); 
      for (i=0; i<nle; ++i){
	fe.coorMap(x,y,z,refFE.xi(i+lcount),refFE.eta(i+lcount),refFE.zeta(i+lcount));
	index = mesh.localEdgeId(ie+1,i+1) - 1;
	supp_x[index]=x;
	supp_y[index]=y;
	supp_z[index]=z;
     }
    }
    gcount+=ne;
    lcount+=nle;
  }
  // Face  Based Dof 
  if (nldpf >0){
    for (ie=0; ie<nV; ++ie){
      fe.updateJac(mesh.volumeList(ie+1)); 
      for (i=0; i<nlf; ++i){
          fe.coorMap(x,y,z,refFE.xi(i+lcount),refFE.eta(i+lcount),refFE.zeta(i+lcount));
          index = mesh.localFaceId(ie+1,i+1)+gcount - 1;
          supp_x[index]=x;
	  supp_y[index]=y;
	  supp_z[index]=z;
      }
    }
    gcount+=nf;
    lcount+=nlf;
  }
  // Volume  Based Dof
  if (nldpV >0 ){
    for (ie=0; ie<nV; ++ie){
      fe.updateJac(mesh.volumeList(ie+1)); 
      fe.coorMap(x,y,z,refFE.xi(lcount),refFE.eta(lcount),refFE.zeta(lcount));
      index = ie + gcount - 1;
      supp_x[index]=x;
      supp_y[index]=y;
      supp_z[index]=z;
    }
  }
  for(i=0;i<num_points_supp;++i){
    ofile << supp_x[i] << " " <<  supp_y[i] << " " <<  supp_z[i] << endl; 
  }

  // connectivity
  // cells_size = nldof*(nV+1);
  cells_size = (nldof+1)*nV;
  
  ofile << endl;
  
  ofile << "CELLS " << nV << " " << cells_size << endl;
  for (i=0;i<nV;++i){
    ofile << nldof << " ";
    for (j=0;j<nldof;++j)   
      ofile << dof.localToGlobal(i+1,j+1)-1 << " ";//damned (C vs) Fortran
    
    ofile << endl;
  }
  
  ofile << endl;
  
  // elements type
  ofile << "CELL_TYPES " << nV << endl;
  for (i=0;i<nV;++i)
    ofile << cell_type << endl;
  ofile << endl;
   ofile << "POINT_DATA " << num_points << endl;
}

/* ! {\tt  wr_vtk_ascii_scalar} is a subrotuine for writing SCALARS unknown in ASCII for VTK
   the file is appended to an existing one (at least you should call before {\tt wr_vtk_ascii_header})
   The {\tt lookup table} is the default one if not differently specified
   In the latter case, it is supposed that in a file (specified by the string look_up_table)
   there is the user defined table.
*/

void wr_vtk_ascii_scalar(string fname, string name, Real* U, int Usize,
			 string look_up_table="default");

void wr_vtk_ascii_vector(string fname, string name, Real* U, int Usize);






//----------------------------------------------------------------------
// obsolete ?
void wr_vtk_ascii_scalar(string fname, string name, vector<Real> U, string look_up_table="default");


/* ! This subroutine considers the "iso" cases, whenever the points of the mesh actually coincide with the ones
 adopted for the computation. In this case, the coordinates of the nodes are ALL immediately available
 and there is no need of recomputing or rereading the other coordinates
 (which on the base of the McKoy algorithm are not explicitly computed so far (see fields...update))
*/
template<typename RegionMesh, typename Dof>
void wr_vtk_ascii_header(string fname, string title, const RegionMesh& mesh, const Dof& dof, UInt cell_type)
{

 ofstream ofile(fname.c_str());
 UInt cells_size; 
 UInt i,j;

  UInt nV=mesh.numVolumes();
  UInt nv=mesh.numVertices();
  UInt nldof=mesh.numLocalVertices();

 ASSERT(ofile,"Error: Output file cannot be open"); 

 ofile << "# vtk DataFile Version " << VTK_VERSION << endl;
 ofile << title << endl;
 ofile << "ASCII" << endl; // capire come scrivere in binario in C++
 ofile << endl;
 ofile << "DATASET UNSTRUCTURED_GRID" << endl;
 ofile << "POINTS " << nv  << " " << "float" << endl; // forse si puo' fare una RTTI sul dato contenuto nella point list
 
 // nodes coordinates (they are all available in the Point List)
 for(i=0;i<nv;++i){
  ofile << mesh.pointList[i].x() << " " <<  mesh.pointList[i].y() << " " <<  mesh.pointList[i].z() << endl; 
 }

 // connectivity
 // cells_size = nldof*(nV+1);
 cells_size = (nldof+1)*nV;
 ofile << endl;

 ofile << "CELLS " << nV << " " << cells_size << endl;
 for (i=0;i<nV;++i){
  ofile << nldof << " ";
  for (j=0;j<nldof;++j)   
    ofile << dof.localToGlobal(i+1,j+1)-1 << " "; //damned (C vs) Fortran
  ofile << endl;
 }
 ofile << endl;
 // elements type
 ofile << "CELL_TYPES " << nV << endl;
 for (i=0;i<nV;++i)
  ofile << cell_type << endl;

 ofile << endl;

}


#endif

