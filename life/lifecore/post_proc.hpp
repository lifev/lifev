/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
/*!
  \file post_proc.h
  \author A. Veneziani
  \date 06/2003
  \version 1.0

  \brief File containing a class with all the methods for a post processing of the solution.
  The basic approach resorts to weak (residual-based) computing of the stresses.

*/
#ifndef _POST_PROC_H
#define _POST_PROC_H
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "lifeV.hpp"
#include "vecUnknown.hpp"

namespace LifeV
{
template<typename Mesh>
class PostProc{

 public:
  PostProc<Mesh>(Mesh& mesh,  CurrentBdFE& feBd, const Dof& dof, UInt ncomp);

  void set_area(CurrentBdFE& feBd, const Dof& dof);
  void set_normal(CurrentBdFE& feBd, const Dof& dof);
  void set_phi(CurrentBdFE& feBd, const Dof& dof);

  void show_bdLtoG();

  void show_area();

  void show_normal();

  void show_phi();

  std::vector< ID > fBdToIn(){return _fBdToIn;}

  Vector compute_sstress(Vector r, UInt ncomp);


 private:
  UInt _nBdDof;
  std::vector<Real> _area; // vector whose ith-component is the area of the patch of the ith-node on the boundary
  std::vector<Real> _normal; // vector with the components of the average  normal vector of each boundary node
  std::vector<Real> _phi;  // vector with \int_{Patch(i)} \phi_i dx where \phi_i is the basis function.

  std::vector< SimpleVect<ID> > _bdLtoG; // for each boundary face, it contains the numbering of the dof of the face
  std::vector< ID > _fBdToIn; // it converts from a local numeration over the boundary faces on the global numeration of the mesh
  Mesh& _mesh;

};

//
// IMPLEMENTATIONS
//
// Construction of _bdLtoG
template <typename Mesh>
PostProc<Mesh>::PostProc(Mesh& mesh, CurrentBdFE& feBd, const Dof& dof, UInt ncomp=1):_mesh(mesh)
{

  typedef  typename Mesh::VolumeShape GeoShape;

  // Some useful local variables, to save some typing
  UInt nDofpV = dof.fe.nbDofPerVertex; // number of Dof per vertices
  UInt nDofpE = dof.fe.nbDofPerEdge;   // number of Dof per edges
  UInt nDofpF = dof.fe.nbDofPerFace;   // number of Dof per faces

  UInt bdnF  = _mesh.numBFaces();    // number of faces on boundary

  typedef typename GeoShape::GeoBShape  GeoBShape;

  UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
  UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges

  UInt nElemV = GeoShape::numVertices; // Number of element's vertices
  UInt nElemE = GeoShape::numEdges;    // Number of element's edges

  UInt nDofFV = nDofpV * nFaceV; // number of vertex's Dof on a face
  UInt nDofFE = nDofpE * nFaceE; // number of edge's Dof on a face

  UInt nDofF  = nDofFV+nDofFE+nDofpF; // number of total Dof on a face

  UInt nDofElemV = nElemV*nDofpV; // number of vertex's Dof on a Element
  UInt nDofElemE = nElemE*nDofpE; // number of edge's Dof on a Element

  SimpleVect<ID> bdltg(nDofF);



  UInt iElAd, iVeEl, iFaEl, iEdEl;
  ID lDof, gDof,auxDof, bDof=1;
  std::vector<ID>::iterator vidit;

  // ===================================================
  // Loop on boundary faces
  // ===================================================
  for (ID ibF=1 ; ibF<=bdnF; ++ibF) {

    iElAd = mesh.boundaryFace(ibF).ad_first();  // id of the element adjacent to the face
    iFaEl = mesh.boundaryFace(ibF).pos_first(); // local id of the face in its adjacent element
    feBd.updateMeas( mesh.boundaryFace(ibF) );  // updating finite element information

    // ===================================================
    // Vertex based Dof
    // ===================================================
    if ( nDofpV ) {

      // loop on face vertices
      for (ID iVeFa=1; iVeFa<=nFaceV; ++iVeFa){
      	iVeEl = GeoShape::fToP(iFaEl,iVeFa); // local vertex number (in element)

	// Loop number of Dof per vertex
	for (ID l=1; l<=nDofpV; ++l) {
	  lDof =   (iVeFa-1) * nDofpV + l ; // local Dof
	  gDof =  dof.localToGlobal( iElAd, (iVeEl-1)*nDofpV + l); // global Dof
          vidit = find(_fBdToIn.begin(),_fBdToIn.end(),gDof);
          if (vidit==_fBdToIn.end()){ // the gDof has been encountered for the first time
            bdltg( lDof ) =  bDof;
            _fBdToIn.push_back(gDof); // local to boundary global on this face
            bDof++;
         }
         else { // the gDof has been already inserted in the _fBdToIn vector
	  auxDof = (ID)((vidit-_fBdToIn.begin()))+1;
	  bdltg( lDof ) =  auxDof; // local to boundary global on this face
	}
       }
      }
    }
    // ===================================================
    // Edge based Dof
    // ===================================================
    if (nDofpE) {
      // loop on face edges
      for (ID iEdFa=1; iEdFa<=nFaceV; ++iEdFa) {
       	iEdEl  = GeoShape::fToE(iFaEl,iEdFa).first; // local edge number (in element)
 	// Loop number of Dof per edge
	for (ID l=1; l<=nDofpE; ++l) {
	  lDof =  nDofFV + (iEdFa-1) * nDofpE + l ; // local Dof
	  gDof =  dof.localToGlobal( iElAd, nDofElemV + (iEdEl-1)*nDofpE + l); // global Dof
          vidit = find(_fBdToIn.begin(),_fBdToIn.end(),gDof);
          if (vidit==_fBdToIn.end()){ // the gDof has been encountered for the first time
            bdltg( lDof ) =  bDof;
            _fBdToIn.push_back(gDof); // local to boundary global on this face
            bDof++;
         }
         else { // the gDof has been already inserted in the _fBdToIn vector
	  auxDof = (ID)(vidit-_fBdToIn.begin())+1;
	  bdltg( lDof ) =  auxDof; // local to boundary global on this face
	}
       }
      }
    }
    // ===================================================
    // Face based Dof
    // ===================================================
	// Loop on number of Dof per face
	for (ID l=1; l<=nDofpF; ++l) {
	  lDof = nDofFE + nDofFV + l; // local Dof
	  gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + (iFaEl-1)*nDofpF + l); // global Dof
          vidit = find(_fBdToIn.begin(),_fBdToIn.end(),gDof);
          if (vidit==_fBdToIn.end()){ // the gDof has been encountered for the first time
            bdltg( lDof ) =  bDof;
            _fBdToIn.push_back(gDof); // local to boundary global on this face
            bDof++;
         }
         else { // the gDof has been already inserted in the _fBdToIn vector
	  auxDof = (ID)(vidit-_fBdToIn.begin())+1;
	  bdltg( lDof ) =  auxDof; // local to boundary global on this face
	}
       }

     _bdLtoG.push_back(bdltg);
  }

  _nBdDof=_fBdToIn.size();
  _area.resize(_nBdDof);
  for (std::vector<Real>::iterator it=_area.begin();it<_area.end();it++)
    *it = 0.0;

  _normal.resize(_nBdDof*NDIM);
  for (std::vector<Real>::iterator it=_normal.begin();it<_normal.end();it++)
    *it = 0.0;

  _phi.resize(_nBdDof);
  for (std::vector<Real>::iterator it=_phi.begin();it<_phi.end();it++)
    *it = 0.0;


}

///////////////////////////////////////////////


// Area of patches on the boundary
template<typename Mesh>
void PostProc<Mesh>::set_area(CurrentBdFE& feBd, const Dof& dof){


  typedef  typename Mesh::VolumeShape GeoShape;

  UInt nDofpV = dof.fe.nbDofPerVertex; // number of Dof per vertices
  UInt nDofpE = dof.fe.nbDofPerEdge;   // number of Dof per edges
  UInt nDofpF = dof.fe.nbDofPerFace;   // number of Dof per faces

  UInt bdnF  = _mesh.numBFaces();    // number of faces on boundary

  typedef typename GeoShape::GeoBShape  GeoBShape;

  UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
  UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges

  UInt nDofFV = nDofpV * nFaceV; // number of vertex's Dof on a face
  UInt nDofFE = nDofpE * nFaceE; // number of edge's Dof on a face

  UInt nDofF  = nDofFV+nDofFE+nDofpF; // number of total Dof on a face

  Real loc_area;

  ID idGlobalDof;


  // ===================================================
  // Loop on boundary faces
  // ===================================================
  for (ID ibF=1 ; ibF<=bdnF; ++ibF) {

    //    iElAd = _mesh.boundaryFace(ibF).ad_first();  // id of the element adjacent to the face
    //    iFaEl = _mesh.boundaryFace(ibF).pos_first(); // local id of the face in its adjacent element
    feBd.updateMeas( _mesh.boundaryFace(ibF) );  // updating finite element information

    loc_area=feBd.measure();
  // Loop on the total Dof per Face
   for (ID idofF=1; idofF<= nDofF; ++idofF) {
     // global dof
     idGlobalDof = _bdLtoG[(UInt)ibF-1][idofF-1]; //
     _area[idGlobalDof-1]+=loc_area;
    }
  }
}


template<typename Mesh>
void PostProc<Mesh>::show_area(){

  ID count=1;

  for (std::vector<Real>::iterator it=_area.begin();it<_area.end();it++){
    std::cout << "Boundary Dof: " << count << ", corresponding to Global Dof: " << _fBdToIn[count-1] << " has patch area: " << *it << std::endl;
    count++;
  }
}

template<typename Mesh>
void PostProc<Mesh>::show_bdLtoG(){

  int count=0;
  std::cout << "***** Post Proc: Bd Local To Global *****" << std::endl;
  std::cout << _bdLtoG.size() << std::endl;
  for (std::vector<SimpleVect<ID> >::iterator it1=_bdLtoG.begin();it1<_bdLtoG.end();it1++) {
    count++;
    std::cout << "Bd Face " << count << std::endl;
    for (SimpleVect<ID>::iterator it2=it1->begin();it2<it1->end();it2++){
      std::cout << *it2 << ",";
   }
    std::cout << std::endl;
  }

 std::cout << "***** Post Proc: From Boundary Faces to Global Dof *****" << std::endl;
  std::cout << _fBdToIn.size() << std::endl;

  for (std::vector<ID>::iterator it3=_fBdToIn.begin();it3<_fBdToIn.end();it3++) {
    std::cout << "Index :" << it3-_fBdToIn.begin() << ", Global Dof: " << *it3 << std::endl;
  }
}

/////////////////////////////////////////////////

///////////////////////////////////////////////


// Normal vectors of patches on the boundary
template<typename Mesh>
void PostProc<Mesh>::set_normal(CurrentBdFE& feBd, const Dof& dof){


  typedef  typename Mesh::VolumeShape GeoShape;

  UInt nDofpV = dof.fe.nbDofPerVertex; // number of Dof per vertices
  UInt nDofpE = dof.fe.nbDofPerEdge;   // number of Dof per edges
  UInt nDofpF = dof.fe.nbDofPerFace;   // number of Dof per faces

  UInt bdnF  = _mesh.numBFaces();    // number of faces on boundary

  typedef typename GeoShape::GeoBShape  GeoBShape;

  UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
  UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges

  UInt nDofFV = nDofpV * nFaceV; // number of vertex's Dof on a face
  UInt nDofFE = nDofpE * nFaceE; // number of edge's Dof on a face

  UInt nDofF  = nDofFV+nDofFE+nDofpF; // number of total Dof on a face

  Real sum;

  ID idGlobalDof;


  // ===================================================
  // Loop on boundary faces
  // ===================================================
  for (ID ibF=1 ; ibF<=bdnF; ++ibF) {

    //    iElAd = _mesh.boundaryFace(ibF).ad_first();  // id of the element adjacent to the face
    //    iFaEl = _mesh.boundaryFace(ibF).pos_first(); // local id of the face in its adjacent element
    feBd.updateMeasNormal( _mesh.boundaryFace(ibF) );  // updating finite element information

    // Loop on the components
    for (int icomp=0; icomp<NDIM;icomp++) {
      sum=0.;
    // Loop on the quadrature points
       for(int l=0; l<feBd.nbQuadPt; ++l) {
	    sum +=  feBd.normal(icomp,l) * feBd.weightMeas(l);
	  }
       for (ID idofF=1; idofF<= nDofF; ++idofF) {
            // global dof
         idGlobalDof = _bdLtoG[(UInt)ibF-1][idofF-1]; //
         _normal[icomp*_nBdDof+idGlobalDof-1]+=sum;
        }
    }
  }
  // Normalization of the averaged normals with the patch area
  for (UInt inorm=0;inorm<_nBdDof;inorm++) {
    Real loc_area=_area[inorm];
    for (int icc=0;icc<NDIM;icc++)
     _normal[icc*_nBdDof+inorm]=_normal[icc*_nBdDof+inorm]/loc_area;

  }


}


template<typename Mesh>
void PostProc<Mesh>::show_normal(){

  ID count=1;

  for (std::vector<Real>::iterator it=_area.begin();it<_area.end();it++){
    std::cout << "Boundary Dof: " << count << ", corresponding to Global Dof: " << _fBdToIn[count-1] << " has patch area: " << *it << std::endl;
    std::cout << "and normal components " ;
    for (int icomp=0; icomp<NDIM; icomp++)
     std::cout << _normal[icomp*_nBdDof+count-1] << " ";

    std::cout << std::endl;
    count++;
  }

  std::cout << "End SHOW NORMAL" << std::endl;
}

//////////////////////////////////////////////////
//
///////////////////////////////////////////////////

// Vector with the integral of the shape functions on the patches on the boundary
template<typename Mesh>
void PostProc<Mesh>::set_phi(CurrentBdFE& feBd, const Dof& dof){


  typedef  typename Mesh::VolumeShape GeoShape;

  UInt nDofpV = dof.fe.nbDofPerVertex; // number of Dof per vertices
  UInt nDofpE = dof.fe.nbDofPerEdge;   // number of Dof per edges
  UInt nDofpF = dof.fe.nbDofPerFace;   // number of Dof per faces

  UInt bdnF  = _mesh.numBFaces();    // number of faces on boundary

  typedef typename GeoShape::GeoBShape  GeoBShape;

  UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
  UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges

  UInt nDofFV = nDofpV * nFaceV; // number of vertex's Dof on a face
  UInt nDofFE = nDofpE * nFaceE; // number of edge's Dof on a face

  UInt nDofF  = nDofFV+nDofFE+nDofpF; // number of total Dof on a face

  Real sum;

  ID idGlobalDof;


  // ===================================================
  // Loop on boundary faces
  // ===================================================
  for (ID ibF=1 ; ibF<=bdnF; ++ibF) {

    feBd.updateMeas( _mesh.boundaryFace(ibF) );  // updating finite element information

     for (ID idofF=1; idofF<= nDofF; ++idofF) {
       sum=0.0;
            // global dof
         idGlobalDof = _bdLtoG[(UInt)ibF-1][idofF-1]; //

     // Loop on the quadrature points
       for(int l=0; l<feBd.nbQuadPt; ++l) {
	    sum +=  feBd.phi((int)(idofF-1),l) * feBd.weightMeas(l);
       }
         _phi[idGlobalDof-1]+=sum;

     }
   }

}


template<typename Mesh>
void PostProc<Mesh>::show_phi(){

  ID count=0;

  for (std::vector<Real>::iterator it=_area.begin();it<_area.end();it++){
    std::cout << "Boundary Dof: " << count+1 << ", corresponding to Global Dof: " << _fBdToIn[count] << " has patch area: " << *it << std::endl;
    std::cout << "and average phi  " << _phi[count]<<std::endl ;
    count++;
  }

  std::cout << "End SHOW PHI" << std::endl;
}



////////////////////////////////
////////////////////////////////
///////////////////////////////
template<typename Mesh>
Vector  PostProc<Mesh>::compute_sstress(Vector r, UInt ncomp){


  ASSERT(ncomp=NDIM,"Error: Shear stress computation possible only for vector unknowns");
  Vector stress(_nBdDof*NDIM);
  Vector nstress(_nBdDof*NDIM);
  Vector sstress(_nBdDof*NDIM);
  ID count=0;
  ID glodof;
  UInt dim=r.size()/ncomp;

  //// (Average) stress computation: from residual to stress
  for (std::vector<Real>::iterator it=_phi.begin();it<_phi.end();it++){
    glodof =  _fBdToIn[count];
    for (UInt ind_comp=0;ind_comp<ncomp;ind_comp++)
      stress[count+ind_comp*_nBdDof]=r[glodof-1+ind_comp*dim]/(*it);//damned conventions trouble : 0 or 1

    count++;
   }

   count=0;
   Real sn=0.;
  ///// Normal stress
  for (count=0;count<_nBdDof;count++){
   sn=0.;
   for (UInt ind_comp=0;ind_comp<ncomp;ind_comp++)
    sn+=stress[count+ind_comp*_nBdDof]*_normal[count+ind_comp*_nBdDof];

   for (UInt ind_comp=0;ind_comp<ncomp;ind_comp++)
    nstress[count+ind_comp*_nBdDof]=sn*_normal[count+ind_comp*_nBdDof];
  }

  // Shear Stress
  sstress = stress - nstress;
  // count=0;
  // for (std::vector<Real>::iterator it=sstress.begin();it<sstress.end();it++)
  // {
  //  *it = stress[count]-nstress[count];
  //  count++;}

  return sstress;
}
}
#endif

