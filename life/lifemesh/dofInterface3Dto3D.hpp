/*!
  \file dofInterface3Dto3D.h
  \brief Class for interfacing dofs between two 3D meshes 
  \version 1.0
  \author M.A. Fernandez
  \date 11/2002

  This file contains the class which may be used to update and hold the connections between the dof
  on two matching meshes.  
*/
#ifndef _DOFINTERFACE3DTO3D_HH
#define _DOFINTERFACE3DTO3D_HH

#include "dofInterfaceBase.hpp"


#include "refFE.hpp"
#include "dof.hpp"
#include <iostream>
#include <map>
#include "tab.hpp"
#include "markers.hpp"
#include "currentBdFE.hpp"
#include <ext/slist> 
using namespace __gnu_cxx;

/*! 
  \class DofInterface3Dto3D

  Base class which holds the conections of the dof in two matching meshes
  
  In order to hold the interface conections the user must give the
  RefFE elements and Dof used in both meshes.  The connections may be
  built by calling the update method. An interpolate method has been
  provided in order to interpolate data at the interface.Finally
  method getInterfaceDof gives the connections.
  
*/
class DofInterface3Dto3D:
  public DofInterfaceBase
{
 public:

  //! Constructor for interfacing Dof of the same type (RefFE)
   /*!
    \param refFe the reference FE used in both meshes
    \param dof1 the Dof object of the mesh in which we want to make the computations
    \param dof2 the Dof object of the mesh which provides de data at the interface
   */
  DofInterface3Dto3D(const RefFE& refFE, const Dof& dof1, const Dof& dof2);
  
  //! Constructor for interfacing Dof of diferent type (RefFE)
  /*!
    \param refFe1 the reference FE used in the mesh in which we want to make the computations
    \param dof1 the Dof object of the mesh in which we want to make the computations
    \param refFe2 the reference FE used in the mesh which provides de data at the interface
    \param dof2 the Dof object of the mesh which provides de data at the interface
   */
  DofInterface3Dto3D(const RefFE& refFE1, const Dof& dof1, const RefFE refFE2, const Dof& dof2);

  //! This method builds the Dof connections at the interface
  /*!
    \param mesh1 the mesh in which we want to make the computations
    \param flag1 the marker of the interface in the mesh1
    \param mesh2 the mesh which provides de data at the interface
    \param flag2 the marker of the interface in the mesh2
    \param tol tolerance for connecting points of both meshes at the interface 
   */
  template<typename Mesh> 
    void update(Mesh& mesh1, const EntityFlag& flag1, Mesh& mesh2, const EntityFlag& flag2, const Real& tol);


  //! This method interpolate data when using different FE
  /*!
    \param mesh2 the mesh which provides de data at the interface
    \param v the data vector on mesh2 holding dofs of type refFE1 
    \param vI the interpolated data vector on mesh2 holding dofs of type refFE2 
    
    \note We should use this method ONLY when the accuracy of the data is less that
          the accuracy of the unknowns on mesh1, i.e., when refFE1 is more accurate
	  than refFE2
   */
  template<typename Mesh, typename VecUnknown> 
    void interpolate(Mesh& mesh2, const VecUnknown& v, VecUnknown& vI);

  //! This method returns the corrresponding dof number of the mesh2 at the interface 
  //! for a specific dof number at the interface in mesh1
  /*!
    \param i a dof number in mesh1 
  */
  // ID getInterfaceDof(const ID& i) const;

 private:

  //! RefFE object used in the mesh in which we want to make the computations
  const RefFE& _refFE1;  

  //! Dof object of the mesh in which we want to make the computations
  const Dof&   _dof1;

  //! RefFE objet used in the mesh which provides de data at the interface
  const RefFE& _refFE2;  

  //! Dof object of the mesh which provides de data at the interface
  const Dof&   _dof2;

  //! Auxiliary Dof object of the mesh which provides de local to global table
  //! when interpolation is used
  Dof _dof;
  
  //! STL list which holds the connections between faces at the interface 
  slist< pair<ID,ID> > _elc;

  //!  Auxiliary STL list which holds the connections between Dof at the interface   
  //! Empty after calling update
  slist< pair<ID,ID> > _locDof;

  //!  STL iterator type for the lists
  typedef slist< pair<ID,ID> >::iterator Iterator;
  
  //! This method builds the connections between faces at the interface (_elc container)
  /*!
    \param mesh1 the mesh in which we want to make the computations
    \param flag1 the marker of the interface in the mesh1
    \param mesh2 the mesh which provides de data at the interface
    \param flag2 the marker of the interface in the mesh2
    \param tol tolerance for connecting points of both meshes at the interface 
   */
  template<typename Mesh>
    void _updateFaceConnections(const Mesh& mesh1, const EntityFlag& flag1, 
				const Mesh& mesh2, const EntityFlag& flag2, const Real& tol);
  
  //! This method builds the connections between Dof at the interface (_locDof container)
  /*!
    \param mesh1 the mesh in which we want to make the computations
    \param dof1 the Dof object of the mesh in which we want to make the computations
    \param mesh2 the mesh which provides de data at the interface
    \param dof2 the Dof object of the mesh which provides de data at the interface
    \param tol tolerance for connecting points of both meshes at the interface 
   */
  template<typename Mesh> 
    void _updateDofConnections(const Mesh& mesh1, const Dof& dof1,
			       const Mesh& mesh2, const Dof& dof2, const Real& tol);  
};


//! Returns true if the vectors v1 and v2 are equal with respect to the tolerance tol 
bool coincide(const KN_<Real>& v1, const KN_<Real>& v2, const Real& tol);


//! Returns true if points (x1,y1,z1) and (x2,y2,z2) are equal with respect to the tolerance tol 
bool coincide(const Real& x1, const Real& y1, const Real& z1, const Real& x2, const Real& y2, const Real& z2, const Real& tol);


 //! This method builds the connections between faces at the interface (_elc container)
  /*!
    \param mesh1 the mesh in which we want to make the computations
    \param flag1 the marker of the interface in the mesh1
    \param mesh2 the mesh which provides de data at the interface
    \param flag2 the marker of the interface in the mesh2
    \param tol tolerance for connecting points of both meshes at the interface 
   */
template<typename Mesh> 
void DofInterface3Dto3D::_updateFaceConnections(const Mesh& mesh1, const EntityFlag& flag1, 
					  const Mesh& mesh2, const EntityFlag& flag2, const Real& tol) {

  UInt bdnF1  = mesh1.numBFaces(); // Number of boundary faces in mesh1
  UInt bdnF2  = mesh2.numBFaces(); // Number of boundary faces in mesh2

  EntityFlag marker1, marker2;
  
  typedef  typename Mesh::FaceShape GeoBShape; // Shape of the faces

  UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices

  KN<Real> v1(nFaceV*nDimensions), v2(nDimensions);
  
  // Loop on boundary faces on mesh1
  for (ID ibF1=1; ibF1<=bdnF1; ++ibF1) {

    // The face marker
    marker1 = mesh1.boundaryFace(ibF1).marker();
    
    // Is the face on the interface?
    if (marker1 == flag1) {
  
      // Loop on face vertices 
      for (ID iVeFa=1; iVeFa <= nFaceV; ++iVeFa) {

	// Loop on vertex coordinates
	for (ID j=1; j<=nDimensions; ++j) 
	  v1[ j-1 + (iVeFa-1)*nDimensions ] = mesh1.boundaryFace(ibF1).point(iVeFa).coordinate(j);
      }
      
      // Loop on boundary faces on mesh2
      for (ID ibF2=1; ibF2<=bdnF2; ++ibF2) {

	// The face marker
	marker2 = mesh2.boundaryFace(ibF2).marker();
	
	// Is the face on the interface?
	if (marker2 == flag2) {  
	 
	  UInt vertexOk=0; // Number of matched vertices

	  // Loop on face vertices 
	  for (ID iVeFa=1; iVeFa <= nFaceV; ++iVeFa) {
	    
	    // Loop on vertex coordinates
	    for (ID j=1; j<=nDimensions; ++j) 
	      v2[ j-1 ] = mesh2.boundaryFace(ibF2).point(iVeFa).coordinate(j);

	    // Loop on face vertices on mesh1 
	    for (ID ivefa=1; ivefa <= nFaceV; ++ivefa) {
	      // Do the vertices match?
	      if ( coincide( v1( SubArray(nDimensions,(ivefa-1)*nDimensions) ),v2,tol) ) {
		++vertexOk;
		break; // Stop loop on face vertices 
	      }
	    }
	  }
	  //! Do the faces match?
	  if (vertexOk == nFaceV) {
	    pair<ID,ID> elc(ibF1,ibF2);
	    _elc.push_front( elc );
	    break; // Stop loop on boundary faces on mesh2
	  }
	}
      }
    }
  }  
}


//! This method builds the connections between Dof at the interface (_locDof container)
/*!
  \param mesh1 the mesh in which we want to make the computations
  \param dof1 the Dof object of the mesh in which we want to make the computations
  \param mesh2 the mesh which provides de data at the interface
  \param dof2 the Dof object of the mesh which provides de data at the interface
  \param tol tolerance for connecting points of both meshes at the interface 
*/
template<typename Mesh> 
void DofInterface3Dto3D::_updateDofConnections(const Mesh& mesh1, const Dof& dof1,
					 const Mesh& mesh2, const Dof& dof2, const Real& tol) {

  typedef  typename Mesh::VolumeShape GeoShape;
  typedef  typename Mesh::FaceShape GeoBShape;
  
  UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
  UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges 
  
  UInt nDofpV1 = _refFE1.nbDofPerVertex; // number of Dof per vertices on mesh1
  UInt nDofpE1 = _refFE1.nbDofPerEdge;   // number of Dof per edges on mesh1
  UInt nDofpF1 = _refFE1.nbDofPerFace;   // number of Dof per faces on mesh1

  UInt nDofpV2 = _refFE2.nbDofPerVertex; // number of Dof per vertices on mesh2
  UInt nDofpE2 = _refFE2.nbDofPerEdge;   // number of Dof per edges on mesh2
  UInt nDofpF2 = _refFE2.nbDofPerFace;   // number of Dof per faces on mesh2

  UInt nElemV = GeoShape::numVertices; // Number of element's vertices 
  UInt nElemE = GeoShape::numEdges;    // Number of element's edges

  UInt nDofFV1 = nDofpV1 * nFaceV; // number of vertex's Dof on a face on mesh1
  UInt nDofFE1 = nDofpE1 * nFaceE; // number of edge's Dof on a face on mesh1

  UInt nDofFV2 = nDofpV2 * nFaceV; // number of vertex's Dof on a face on mesh2
  UInt nDofFE2 = nDofpE2 * nFaceE; // number of edge's Dof on a face on mesh2

  UInt nDofElemV1 = nElemV*nDofpV1; // number of vertex's Dof on a Element on mesh1
  UInt nDofElemE1 = nElemE*nDofpE1; // number of edge's Dof on a Element on mesh1

  UInt nDofElemV2 = nElemV*nDofpV2; // number of vertex's Dof on a Element on mesh2
  UInt nDofElemE2 = nElemE*nDofpE2; // number of edge's Dof on a Element on mesh2

  ID iElAd1, iVeEl1, iFaEl1, iEdEl1, iElAd2, iVeEl2, iFaEl2, iEdEl2, lDof1, lDof2, gDof1, gDof2;

  Real x1, x2, y1, y2, z1, z2;

  bool test = false;

  CurrentBdFE feBd1( _refFE1.boundaryFE(), mesh1.getGeoMap().boundaryMap() ); 
  CurrentBdFE feBd2( _refFE2.boundaryFE(), mesh2.getGeoMap().boundaryMap() );
  
  // Loop on faces at the interface (matching faces)
  for (Iterator i=_elc.begin(); i != _elc.end(); ++i) {
    
    feBd1.update( mesh1.boundaryFace(i->first) );  // Updating face information on mesh1
    feBd2.update( mesh2.boundaryFace(i->second) );  // Updating face information on mesh2
    
    iElAd1 = mesh1.boundaryFace(i->first).ad_first();  // id of the element adjacent to the face (mesh1) 
    iElAd2 = mesh2.boundaryFace(i->second).ad_first();  // id of the element adjacent to the face (mesh2)
    
    iFaEl1 = mesh1.boundaryFace(i->first).pos_first(); // local id of the face in its adjacent element (mesh1)
    iFaEl2 = mesh2.boundaryFace(i->second).pos_first(); // local id of the face in its adjacent element (mesh2)
        
    // Vertex based Dof on mesh1
    if ( nDofpV1 ) { 
      
      // loop on face vertices (mesh1)
      for (ID iVeFa1=1; iVeFa1<=nFaceV; ++iVeFa1){
	
	iVeEl1 = GeoShape::fToP(iFaEl1,iVeFa1); // local vertex number (in element)
	
	// Loop number of Dof per vertex (mesh1)
	for (ID l=1; l<=nDofpV1; ++l) {
	  lDof1 = (iVeFa1-1) * nDofpV1 + l ; // local Dof 
	  feBd1.coorMap(x1, y1, z1, feBd1.refFE.xi(lDof1-1), feBd1.refFE.eta(lDof1-1) ); // Nodal coordinates on the current face (mesh1)
	    
	  // loop on face vertices (mesh2)
	  for (ID iVeFa2=1; iVeFa2<=nFaceV; ++iVeFa2){
	    
	    iVeEl2 = GeoShape::fToP(iFaEl2,iVeFa2); // local vertex number (in element)
	    
	    // Loop on number of Dof per vertex (mesh2)
	    for (ID k=1; k<=nDofpV2; ++k) {
	      lDof2 =   (iVeFa2-1) * nDofpV2 + k ; // local Dof 
	      feBd2.coorMap(x2, y2, z2, feBd2.refFE.xi(lDof2-1), feBd2.refFE.eta(lDof2-1) ); // Nodal coordinates on the current face (mesh2)
	      	    
	      // Do the nodal points match?
	      if ( test=coincide(x1,y1,z1,x2,y2,z2,tol) ) {
		gDof1 = dof1.localToGlobal(iElAd1,(iVeEl1-1)*nDofpV1 + l); // Global Dof on mesh1
		gDof2 = dof2.localToGlobal(iElAd2,(iVeEl2-1)*nDofpV2 + k); // Global Dof on mesh2
		pair<ID,ID> locDof(gDof1, gDof2); 
		_locDof.push_front( locDof ); // Updating the list of dof connections
		break;
	      }
	    }
	    // Exit the loop on face vertices on mesh2?
	    if ( test ) {
	      test = false;
	      break;
	    }
	  }    
	}
      }
    }
       
    // Edge based Dof on mesh1
    if (nDofpE1) { 
	
      // loop on face edges (mesh1)
      for (ID iEdFa1=1; iEdFa1<=nFaceE; ++iEdFa1) {
	
	iEdEl1  = GeoShape::fToE(iFaEl1,iEdFa1).first; // local edge number (in element)
	
	// Loop number of Dof per edge (mesh1)
	for (ID l=1; l<=nDofpE1; ++l) {
	  
	  lDof1 =  nDofFV1 + (iEdFa1-1) * nDofpE1 + l ; // local Dof 
	  feBd1.coorMap(x1, y1, z1, feBd1.refFE.xi(lDof1-1), feBd1.refFE.eta(lDof1-1) ); // Nodal coordinates on the current face (mesh1)
	  
	  // loop on face edges (mesh2)
	  for (ID iEdFa2=1; iEdFa2<=nFaceV; ++iEdFa2) {
	    
	    iEdEl2 = GeoShape::fToE(iFaEl2,iEdFa2).first; // local edge number (in element)
	    
	    // Loop number of Dof per edge (mesh1)
	    for (ID k=1; k<=nDofpE2; ++k) {
	      
	      lDof2 = nDofFV2 + (iEdFa2-1) * nDofpE2 + k; // local Dof 
	      feBd2.coorMap(x2, y2, z2, feBd2.refFE.xi(lDof2-1), feBd2.refFE.eta(lDof2-1) ); // Nodal coordinates on the current face (mesh2)
 	    
	      // Do the nodal points match?
	      if ( test=coincide(x1,y1,z1,x2,y2,z2,tol) ) {
		gDof1 = dof1.localToGlobal(iElAd1, nDofElemV1 + (iEdEl1-1)*nDofpE1 + l); // Global Dof on mesh1
		gDof2 = dof2.localToGlobal(iElAd2, nDofElemV2 + (iEdEl2-1)*nDofpE2 + k ); // Global Dof on mesh2
		pair<ID,ID> locDof(gDof1, gDof2); 
		_locDof.push_front( locDof ); // Updating the list of dof connections
		break;
		}
	    }
	    // Exit the loop on face edges on mesh2?
	    if ( test ) {
	      test = false;
	      break;
	    }
	  }
	}
      }
    }
        
    // Face based Dof on mesh1
    for (ID l=1; l<=nDofpF1; ++l) {  
      lDof1 = nDofFE1 + nDofFV1 + l; // local Dof 
      feBd1.coorMap(x1, y1, z1, feBd1.refFE.xi(lDof1-1), feBd1.refFE.eta(lDof1-1) ); // Nodal coordinates on the current face (mesh1)
      
      for (ID k=1; k<=nDofpF2; ++k) {  
	lDof2 = nDofFE2 + nDofFV2 + k; // local Dof 
	feBd2.coorMap(x2, y2, z2, feBd2.refFE.xi(lDof2-1), feBd2.refFE.eta(lDof2-1) ); // Nodal coordinates on the current face (mesh2)

	// Do the nodal points match?
	if ( coincide(x1,y1,z1,x2,y2,z2,tol) ) {
	  gDof1 = dof1.localToGlobal(iElAd1, nDofElemE1 + nDofElemV1 + (iFaEl1-1)*nDofpF1 + l); // Global Dof in mesh1 
	  gDof2 = dof2.localToGlobal(iElAd2, nDofElemE2 + nDofElemV2 + (iFaEl2-1)*nDofpF2 + k); // Global Dof in mesh2
	  pair<ID,ID> locDof(gDof1, gDof2);
	  _locDof.push_front( locDof ); // Updating the list of dof connections
	  break;
	}
      }   
    }
  } 

  // Updating the map containter with the connections
  for (Iterator i=_locDof.begin(); i!=_locDof.end(); ++i) {
    _locDofMap[i->first] = i->second;
  }
  
  // Saving memory
  _locDof.clear();
}


//! This method builds the Dof connections at the interface
/*!
  \param mesh1 the mesh in which we want to make the computations
  \param flag1 the marker of the interface in the mesh1
  \param mesh2 the mesh which provides de data at the interface
  \param flag2 the marker of the interface in the mesh2
  \param tol tolerance for connecting points of both meshes at the interface 
*/
template<typename Mesh>  
void DofInterface3Dto3D::update(Mesh& mesh1, const EntityFlag& flag1, 
			  Mesh& mesh2, const EntityFlag& flag2, const Real& tol) {

  // Updating face connections at the interface
  _updateFaceConnections(mesh1,flag1,mesh2,flag2,tol);
  
  if ( _refFE1.nbDof > _refFE2.nbDof ) {
    // Update of the Dof connections when we need interpolation
    _dof.update(mesh2); // Building auxiliary dof  
    _updateDofConnections(mesh1, _dof1, mesh2, _dof, tol); // Update of the Dof connections
  }
  else  
    // Update of the Dof connections without inperpolation
    _updateDofConnections(mesh1, _dof1, mesh2, _dof2, tol);
}


//! This method interpolate data when using different FE
/*!
  \param mesh2 the mesh which provides de data at the interface
  \param v the data vector on mesh2 holding dofs of type refFE1 
  \param vI the interpolated data vector on mesh2 holding dofs of type refFE2 
  
  \note We should use this method ONLY when the acuracy of the data is less that
  the accuracy of the unknowns on mesh1, i.e., when refFE1 is more accurate
  than refFE2
*/
template<typename Mesh, typename VecUnknown> 
void DofInterface3Dto3D::interpolate(Mesh& mesh2, const VecUnknown& v, VecUnknown& vI) {
  
  typedef  typename Mesh::VolumeShape GeoShape; // Element shape
  typedef  typename Mesh::FaceShape GeoBShape;  // Face Shape
  
  UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
  UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges 
  
  UInt nDofpV = _refFE1.nbDofPerVertex; // number of Dof per vertices
  UInt nDofpE = _refFE1.nbDofPerEdge;   // number of Dof per edges
  UInt nDofpF = _refFE1.nbDofPerFace;   // number of Dof per faces
  
  UInt nElemV = GeoShape::numVertices; // Number of element's vertices 
  UInt nElemE = GeoShape::numEdges;    // Number of element's edges
  
  UInt nDofElem = _refFE2.nbDof; // Number of Dof per element in the lowDof mesh
  
  UInt nDofElemV = nElemV*nDofpV; // number of vertex's Dof on a Element
  UInt nDofElemE = nElemE*nDofpE; // number of edge's Dof on a Element
  
  ID ibF, iElAd, iFaEl, iVeEl, lDof, iEdEl;
  
  ID nbComp = v.nbcomp(); // Number of components of the data vector

  Real x, y, z, sum;

  KN<Real> vLoc( nDofElem*nbComp );
 
  
  // Loop on faces at the interface (matching faces)
  for (Iterator i=_elc.begin(); i != _elc.end(); ++i) {
    
    ibF = i->second; // Face number at the interface
        
    iElAd = mesh2.boundaryFace(ibF).ad_first();  // id of the element adjacent to the face 
    iFaEl = mesh2.boundaryFace(ibF).pos_first(); // local id of the face in its adjacent element
 
    // Updating the local dof of the data vector in the adjacent element
    for (UInt icmp=0; icmp < nbComp; ++icmp)
      for (ID idof=0; idof < nDofElem; ++idof) 
	vLoc( icmp*nDofElem + idof ) = (v.vec())( icmp*_dof2.numTotalDof() + _dof2.localToGlobal(iElAd,idof+1)-1 );

    // Vertex based Dof 
    if ( nDofpV ) { 
      
      // loop on face vertices 
      for (ID iVeFa=1; iVeFa<=nFaceV; ++iVeFa){
	
	iVeEl = GeoShape::fToP(iFaEl,iVeFa); // local vertex number (in element)
	
	// Loop number of Dof per vertex
	for (ID l=1; l<=nDofpV; ++l) {
	  lDof = (iVeEl-1)*nDofpV + l; // Local dof in the adjacent Element

	  // Nodal coordinates
	  x = _refFE1.xi(lDof-1);
	  y = _refFE1.eta(lDof-1);
	  z = _refFE1.zeta(lDof-1);
	  
	  // Loop on data vector components
	  for (UInt icmp=0; icmp < nbComp; ++icmp) {

	    // Interpolating data at the nodal point
	    sum = 0;
	    for (ID idof=0; idof < nDofElem; ++idof) // Loop on local Dof on the adjacent element 
	      sum += vLoc(icmp*nDofElem + idof) * _refFE2.phi(idof, x, y, z);
	   
	    // Updating interpolating vector
	    (vI.vec())( icmp*_dof.numTotalDof() + _dof.localToGlobal(iElAd,lDof) - 1 ) = sum; 
	  }
	}
      }
    }
    
    // Edge based Dof 
    if (nDofpE) { 
	
      // loop on face edges 
      for (ID iEdFa=1; iEdFa<=nFaceE; ++iEdFa) {
	
	iEdEl  = GeoShape::fToE(iFaEl,iEdFa).first; // local edge number (in element)
	
	// Loop number of Dof per edge
	for (ID l=1; l<=nDofpE; ++l) {
	  lDof = nDofElemV + (iEdEl-1)*nDofpE + l; // Local dof in the adjacent Element
	 
	  // Nodal coordinates
	  x = _refFE1.xi(lDof-1);
	  y = _refFE1.eta(lDof-1);
	  z = _refFE1.zeta(lDof-1);
	  
	  // Loop on data vector components
	  for (UInt icmp=0; icmp < nbComp; ++icmp) {  

	    // Interpolating data at the nodal point
	    sum = 0;
	    for (ID idof=0; idof < nDofElem; ++idof)  // Loop on local Dof on the adjacent element 
	      sum += vLoc(icmp*nDofElem + idof) * _refFE2.phi(idof, x, y, z);
	  
	    // Updating interpolating vector
	    (vI.vec())( icmp*_dof.numTotalDof() + _dof.localToGlobal(iElAd,lDof) - 1 ) = sum;
	  }   
	}
      }
    }  
    
    // Loop on number of Dof per face
    for (ID l=1; l<=nDofpF; ++l) {  
      lDof = nDofElemE + nDofElemV + (iFaEl-1)*nDofpF + l; // Local dof in the adjacent Element
     
      // Nodal coordinates
      x = _refFE1.xi(lDof-1);
      y = _refFE1.eta(lDof-1);
      z = _refFE1.zeta(lDof-1);
      
      // Loop on data vector components
      for (UInt icmp=0; icmp < nbComp; ++icmp) {
	
	// Interpolating data at the nodal point
	sum = 0;
	for (ID idof=0; idof < nDofElem; ++idof) // Loop on local Dof on the adjacent element 
	  sum += vLoc(icmp*nDofElem + idof) * _refFE2.phi(idof, x, y, z);

	// Updating interpolating vector
	(vI.vec())( icmp*_dof.numTotalDof() + _dof.localToGlobal(iElAd,lDof) - 1) = sum;
      }      
    }
  }
}    
 
#endif
