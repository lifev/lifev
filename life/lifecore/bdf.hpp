/*!
  \file bdf.h
  \author A. Veneziani
  \date 04/2003 
  \version 1.0

  \brief File containing a class for an easy handling of different order time discretizations/extrapolations
         BDF based 

*/
#ifndef _BDF_H
#define _BDF_H
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "lifeV.hpp"
#include "vecUnknown.hpp"

#define MAX_ORDER 3

typedef Real (*Funct)(const Real&, const Real&, const Real&, const Real&, const ID&);

class Bdf
{

 public:
  // !Constructor
  Bdf(const UInt n);
  // ! Inizialize all the entries of the unknown vector to be derived with the vector u0 (duplicated)
  void initialize_unk(Vector u0);
  // ! Inizialize all the entries of the unknown vector to be derived with a set of vectors uv0
  void initialize_unk(vector<Vector > uv0);
  // ! Initialize alle the entries of the unknonwn vectors with a given function
template<typename Mesh,typename RefFE, typename CurrFE, typename Dof> 
  void initialize_unk(const Funct& u0,Mesh& mesh,RefFE& refFE, CurrFE& currFE, Dof& dof, Real t0, Real dt, UInt nbComp);


  // ! Update the vectors of the previous time steps by shifting on the right the old values
  void shift_right(Vector u_curr);

  // ! Carry out the time derivative (BDF of order q)
  Vector time_der(Real dt);

  // ! Carry out the time derivative (BDF of order q) whenever the time step is considered elsewhere
  Vector time_der();

  // ! Carry out the time extrapolation (BDF based)
  Vector extrap();

  // ! Return the i-th coefficient of the time derivative 
  double coeff_der(UInt i); 

  // ! Return the i-th coefficient of the time extrapolation 
  double coeff_ext(UInt i); 

  vector<Vector> unk();

  void showMe();

  ~Bdf();

 private:
  // ! Order of the BDF derivative/extrapolation: the time-derivative coefficients vector has size n+1,
  // ! the extrapolation vector has size n
  UInt _n;
  // ! Size of the unknown vector
  UInt _s;
  // ! Coefficients of the time bdf discretization: u_t (t^{n+1}) = alpha(0) u^{n+1}-alpha(1) u^{n} - alpha(2) u^{n-1} ...
  Vector _alpha; 
  // ! Coefficients of the  extrapolation: u(t^{n+1}) = beta(0) u^{n}+beta(1) u^{n-1} + beta(2) u^{n-2} + ...
  Vector _beta;

  vector<Vector > _unk; // Unknown solutions stored
};


///
// template implementations
//

template<typename Mesh, typename RefFE, typename CurrFE, typename Dof>
void Bdf::initialize_unk(const Funct& u0,Mesh& mesh,RefFE& refFE, CurrFE& currFE, Dof& dof, Real t0, Real dt, UInt nbComp=1)
  // The array of initial conditions needed by the selected BDF is initialized as follows 
  // _unk=[ u0(t0), u0(t0-dt), u0(t0-2*dt), ...]
  // For the space dependence of the initial conditions we need informations on:
  // 1) the mesh (coordinates of points)
  // 2) the reference and the current FE (basis functions)
  // 3) is it a vector or a scalar problem ? bdf doesn't know it
  // 4) which is the initial time (t0) and the time step (for solutions before the initial instant)
  // Based on the NavierStokesHandler::initialize by M. Fernandez 
{
  typedef  typename Mesh::VolumeShape GeoShape; // Element shape
   
  UInt nDofpV  = refFE.nbDofPerVertex; // number of Dof per vertex
  UInt nDofpE  = refFE.nbDofPerEdge;   // number of Dof per edge
  UInt nDofpF  = refFE.nbDofPerFace;   // number of Dof per face
  UInt nDofpEl = refFE.nbDofPerVolume; // number of Dof per Volume
  
  UInt nElemV = GeoShape::numVertices; // Number of element's vertices 
  UInt nElemE = GeoShape::numEdges;    // Number of element's edges
  UInt nElemF = GeoShape::numFaces;    // Number of element's faces
  
  UInt nDofElemV = nElemV*nDofpV; // number of vertex's Dof on a Element
  UInt nDofElemE = nElemE*nDofpE; // number of edge's Dof on a Element
  UInt nDofElemF = nElemF*nDofpF; // number of face's Dof on a Element
    
  vector< Vector >::iterator iter=_unk.begin();
  vector< Vector >::iterator iter_end=_unk.end();


  UInt size_comp=dof.numTotalDof();
  _s=size_comp*nbComp;// Inizialization of the dimension of the vector

  Vector aux(_s);
  aux=0.0;


  for (  iter=_unk.begin() ; iter != iter_end; iter++){
    *iter=aux;
  }

  Real x, y, z;
  UInt backtime; 
  ID lDof;

  // Loop on elements of the mesh
  for (ID iElem=1; iElem <= mesh.numVolumes(); ++iElem){
       
    currFE.updateJac(mesh.volume(iElem) );

    // Vertex based Dof 
    if ( nDofpV ) { 
      
      // loop on element vertices 
      for (ID iVe=1; iVe<=nElemV; ++iVe){
	
	// Loop number of Dof per vertex
	for (ID l=1; l<=nDofpV; ++l) {
	  lDof = (iVe-1)*nDofpV + l; // Local dof in this element
	  
	  // Nodal coordinates
	  currFE.coorMap(x, y, z, currFE.refFE.xi(lDof-1), currFE.refFE.eta(lDof-1), currFE.refFE.zeta(lDof-1));  
	  
	  // Loop on data vector components ****************
          backtime=0;
          for (vector<Vector>::iterator it=_unk.begin();it<_unk.end();it++){
	    for (UInt icmp=0; icmp < nbComp; ++icmp){ 
	      //              cout << "comp " << icmp*size_comp + dof.localToGlobal(iElem,lDof) - 1 << endl;
	     (*it)( icmp*size_comp + dof.localToGlobal(iElem,lDof) - 1 ) = u0(t0-backtime*dt,x,y,z,icmp+1);
             backtime++;
	    }
	  }
	}
      }
    }

    // Edge based Dof 
    if (nDofpE) { 
	
      // loop on element edges 
      for (ID iEd=1; iEd <=nElemE; ++iEd) {
	
	// Loop number of Dof per edge
	for (ID l=1; l<=nDofpE; ++l) {
	  lDof = nDofElemV + (iEd-1)*nDofpE + l; // Local dof in the adjacent Element
	
	  // Nodal coordinates
	  currFE.coorMap(x, y, z, currFE.refFE.xi(lDof-1), currFE.refFE.eta(lDof-1), currFE.refFE.zeta(lDof-1));
	
	  // Loop on data vector components
          backtime=0;
          for (vector<Vector>::iterator it=_unk.begin();it<_unk.end();it++){
	    for (UInt icmp=0; icmp < nbComp; ++icmp){ 
	     (*it)( icmp*size_comp + dof.localToGlobal(iElem,lDof) - 1 ) = u0(t0-backtime*dt,x,y,z,icmp+1);
             backtime++;
	    }
	  }
	}
      }
    }  
    
    // Face based Dof 
    if (nDofpF) { 
      
      // loop on element faces
      for (ID iFa=1; iFa <=nElemF; ++iFa) {
	
	// Loop on number of Dof per face
	for (ID l=1; l<=nDofpF; ++l) {
	  
	  lDof = nDofElemE + nDofElemV + (iFa-1)*nDofpF + l; // Local dof in the adjacent Element
     
	  // Nodal coordinates
	  currFE.coorMap(x, y, z, currFE.refFE.xi(lDof-1), currFE.refFE.eta(lDof-1), currFE.refFE.zeta(lDof-1));
		  
	  // Loop on data vector components
          backtime=0;
          for (vector<Vector>::iterator it=_unk.begin();it<_unk.end();it++){
	    for (UInt icmp=0; icmp < nbComp; ++icmp){ 
	     (*it)( icmp*size_comp + dof.localToGlobal(iElem,lDof) - 1 ) = u0(t0-backtime*dt,x,y,z,icmp+1);
             backtime++;
	    }
	  }
	}
      }
    }

    // Element based Dof 
    // Loop on number of Dof per Element
    for (ID l=1; l<=nDofpEl; ++l) {
      lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element
        
      // Nodal coordinates
      currFE.coorMap(x, y, z, currFE.refFE.xi(lDof-1), currFE.refFE.eta(lDof-1), currFE.refFE.zeta(lDof-1));
	      
      // Loop on data vector components
          backtime=0;
          for (vector<Vector>::iterator it=_unk.begin();it<_unk.end();it++){
	    for (UInt icmp=0; icmp < nbComp; ++icmp){ 
	     (*it)( icmp*size_comp + dof.localToGlobal(iElem,lDof) - 1 ) = u0(t0-backtime*dt,x,y,z,icmp+1);
             backtime++;
	    }
	  }
    }
  }
  return;
}

#endif
