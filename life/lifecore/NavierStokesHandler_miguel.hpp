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
  \file NavierStokesHandler.h
  \author M.A. Fernandez
  \date 06/2003 
  \version 1.0

  \brief This file contains an abstract class for NavierStokes solvers.  

*/

#ifndef _NAVIERSTOKESHANDLER_H_
#define _NAVIERSTOKESHANDLER_H_

#include <cmath>
#include <sstream> 
#include <cstdlib>


#include "dataNavierStokes.hpp"
#include "dataAztec.hpp"
#include "refFE.hpp"
#include "dof.hpp"
#include "lifeV.hpp"
#include "medit_wrtrs.hpp"
#include "bcCond.hpp"


/*! 
  \class NavierStokesHandler

  Abstract class which defines the general structure of a NavierStokes solver.
  For each new NavierStokes solver  we have to implement the corresponding timeAdvance and an iterate methods 
  
*/

template <typename Mesh>
class NavierStokesHandler:
public DataNavierStokes<Mesh> { 
 
 public:

  typedef Real (*Function)(const Real&, const Real&, const Real&, const Real&, const ID&);

  //! Constructor   
  /*!
    \param data_file GetPot data file
    \param refFE_u reference FE for the velocity
    \param refFE_p reference FE for the pressure
    \param Qr_u volumic quadrature rule for the velocity
    \param bdQr_u surface quadrature rule for the velocity 
    \param Qr_p volumic quadrature rule for the pressure
    \param bdQr_p surface quadrature rule for the pressure
    \param BCh_u boundary conditions for the velocity
  */
  NavierStokesHandler(const GetPot& data_file,  const RefFE& refFE_u, 
		const RefFE& refFE_p, const QuadRule& Qr_u, const QuadRule& bdQr_u, 
		const QuadRule& Qr_p, const QuadRule& bdQr_p, BC_Handler& BCh_u);
    
  //! Sets initial condition for the velocity
  void initialize(const Function& u0); 

  //! Update the right  hand side  for time advancing   
  /*! 
    \param source volumic source  
    \param time present time
  */
  virtual void timeAdvance(const Function source, const Real& time) = 0; 

  //! Update convective term, bc treatment and solve the linearized ns system
  virtual void iterate() = 0;

  //! Returns the velocity vector
  PhysVectUnknown<Vector>& u();
  
  //! Returns the pressure
  ScalUnknown<Vector>& p();

  //! Returns the velocity Dof 
  const Dof& uDof() const;

  //! Returns the pressure Dof 
  const Dof& pDof() const;

  //! Postprocessing
  void postProcess();

  //! Do nothing destructor
  virtual ~NavierStokesHandler() {}
  
 protected:
  
  //! Reference FE for the velocity
  const RefFE& _refFE_u;
  
  //! Reference FE for the pressure
  const RefFE& _refFE_p;
 
  //! The Dof object associated with the velocity
  Dof _dof_u;

  //! The Dof object associated with the pressure
  Dof _dof_p;

  //! The number of total velocity dofs  
  UInt _dim_u;

  //! The number of total pressure dofs  
  UInt _dim_p;

  //! Quadrature rule for velocity volumic elementary computations
  const QuadRule& _Qr_u;
  
  //! Quadrature rule for velocity surface elementary computations
  const QuadRule& _bdQr_u;

  //! Quadrature rule for pressure volumic elementary computations
  const QuadRule& _Qr_p;
  
  //! Quadrature rule for pressure surface elementary computations
  const QuadRule& _bdQr_p;

  //! Current FE for the velocity u
  CurrentFE _fe_u;
  CurrentBdFE _feBd_u;

  //! Current FE for the pressure p
  CurrentFE _fe_p;

  //! The velocity
  PhysVectUnknown<Vector> _u;

  //! The pressure
  ScalUnknown<Vector> _p;
   
  //! The BC handler
  BC_Handler& _BCh_u;

  //! The actual time
  Real _time;

  //! Aux. var. for PostProc
  UInt _count;
};



//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh> 
NavierStokesHandler<Mesh>::
NavierStokesHandler(const GetPot& data_file,  const RefFE& refFE_u, 
		    const RefFE& refFE_p, const QuadRule& Qr_u, const QuadRule& bdQr_u, 
		    const QuadRule& Qr_p, const QuadRule& bdQr_p, BC_Handler& BCh_u):
     DataNavierStokes<Mesh>(data_file),    
     _refFE_u(refFE_u),
     _refFE_p(refFE_p),
     _dof_u(_mesh,_refFE_u),
     _dof_p(_mesh,_refFE_p),
     _dim_u(_dof_u.numTotalDof()),
     _dim_p(_dof_p.numTotalDof()), 
     _Qr_u(Qr_u),
     _bdQr_u(bdQr_u),
     _Qr_p(Qr_p),
     _bdQr_p(bdQr_p),
     _fe_u(_refFE_u,_mesh.getGeoMap(),_Qr_u),
     _feBd_u(_refFE_u.boundaryFE(),_mesh.getGeoMap().boundaryMap(),_bdQr_u),
     _fe_p(_refFE_p,_mesh.getGeoMap(),_Qr_p),
     _u(_dim_u),
     _p(_dim_p),
     _BCh_u(BCh_u),
     _time(0),
     _count(0) {}


// Returns the velocity vector
template<typename Mesh> PhysVectUnknown<Vector>& 
NavierStokesHandler<Mesh>::u() {
  return _u;
}

// Returns the pressure
template<typename Mesh> ScalUnknown<Vector>&
NavierStokesHandler<Mesh>::p() {
  return _p;
}

// Returns the velocity Dof 
template<typename Mesh> const Dof& 
NavierStokesHandler<Mesh>::uDof() const {
  return _dof_u;
}

// Returns the pressure Dof 
template<typename Mesh>  const Dof& 
NavierStokesHandler<Mesh>::pDof() const {
  return _dof_p;
}

// Postprocessing 
template<typename Mesh>  void 
NavierStokesHandler<Mesh>::postProcess() {
  ostringstream index;
  string name;
 
  ++_count;
 
  if (fmod(float(_count),float(_verbose)) == 0.0) {
    cout << "  o-  Post-processing \n";
    index << (_count/_verbose);
 
    switch( index.str().size() ) {
    case 1:
      name = "00"+index.str();
      break;
    case 2:
      name = "0"+index.str();
      break;
    case 3:
      name = index.str();
      break;
    }
 
 
    wr_medit_ascii_scalar("press."+name+".bb",_p.giveVec(),_p.size());
    wr_medit_ascii_scalar("vel_x."+name+".bb",_u.giveVec(),_mesh.numVertices());
    wr_medit_ascii_scalar("vel_y."+name+".bb",_u.giveVec() + _dim_u,_mesh.numVertices());
    wr_medit_ascii_scalar("vel_z."+name+".bb",_u.giveVec() + 2*_dim_u,_mesh.numVertices());
    // wr_medit_ascii_vector("veloc."+name+".bb",_u.giveVec(),_mesh.numVertices(),_dim_u);
    system(("ln -s "+_mesh_file+" press."+name+".mesh").data());
    system(("ln -s "+_mesh_file+" vel_x."+name+".mesh").data());
    system(("ln -s "+_mesh_file+" vel_y."+name+".mesh").data());
    system(("ln -s "+_mesh_file+" vel_z."+name+".mesh").data());
    // system(("ln -s "+_mesh_file+" veloc."+name+".mesh").data());
  }
}


// Sets the initial condition
template<typename Mesh> void 
NavierStokesHandler<Mesh>::initialize(const Function& u0) {
  
  // Initialize pressure
  _p.vec()=0.0;

  // Initialize velocity

  typedef  typename Mesh::VolumeShape GeoShape; // Element shape
   
  UInt nDofpV  = _refFE_u.nbDofPerVertex; // number of Dof per vertex
  UInt nDofpE  = _refFE_u.nbDofPerEdge;   // number of Dof per edge
  UInt nDofpF  = _refFE_u.nbDofPerFace;   // number of Dof per face
  UInt nDofpEl = _refFE_u.nbDofPerVolume; // number of Dof per Volume
  
  UInt nElemV = GeoShape::numVertices; // Number of element's vertices 
  UInt nElemE = GeoShape::numEdges;    // Number of element's edges
  UInt nElemF = GeoShape::numFaces;    // Number of element's faces
  
  UInt nDofElemV = nElemV*nDofpV; // number of vertex's Dof on a Element
  UInt nDofElemE = nElemE*nDofpE; // number of edge's Dof on a Element
  UInt nDofElemF = nElemF*nDofpF; // number of face's Dof on a Element
    
  ID nbComp = _u.nbcomp(); // Number of components of the mesh velocity

  Real x, y, z;
 
  ID lDof;

  // Loop on elements of the mesh
  for (ID iElem=1; iElem <= _mesh.numVolumes(); ++iElem) {
       
    _fe_u.updateJac( _mesh.volume(iElem) );

    // Vertex based Dof 
    if ( nDofpV ) { 
      
      // loop on element vertices 
      for (ID iVe=1; iVe<=nElemV; ++iVe){
	
	// Loop number of Dof per vertex
	for (ID l=1; l<=nDofpV; ++l) {
	  lDof = (iVe-1)*nDofpV + l; // Local dof in this element
	  
	  // Nodal coordinates
	  _fe_u.coorMap(x, y, z, _fe_u.refFE.xi(lDof-1), _fe_u.refFE.eta(lDof-1), _fe_u.refFE.zeta(lDof-1));  
	  
	  // Loop on data vector components
	  for (UInt icmp=0; icmp < nbComp; ++icmp)
       	    _u.vec()( icmp*_dim_u + _dof_u.localToGlobal(iElem,lDof) - 1 ) = u0(0.0,x,y,z,icmp+1);  
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
	  _fe_u.coorMap(x, y, z, _fe_u.refFE.xi(lDof-1), _fe_u.refFE.eta(lDof-1), _fe_u.refFE.zeta(lDof-1));
	
	  // Loop on data vector components
	  for (UInt icmp=0; icmp < nbComp; ++icmp) 
	    _u.vec()( icmp*_dim_u + _dof_u.localToGlobal(iElem,lDof) - 1 ) = u0(0.0,x,y,z,icmp+1);
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
	  _fe_u.coorMap(x, y, z, _fe_u.refFE.xi(lDof-1), _fe_u.refFE.eta(lDof-1), _fe_u.refFE.zeta(lDof-1));
		  
	  // Loop on data vector components
	  for (UInt icmp=0; icmp < nbComp; ++icmp) 
	    _u.vec()( icmp*_dim_u + _dof_u.localToGlobal(iElem,lDof) - 1) = u0(0.0,x,y,z,icmp+1);   
	}
      }
    }

    // Element based Dof 
    // Loop on number of Dof per Element
    for (ID l=1; l<=nDofpEl; ++l) {
      lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element
        
      // Nodal coordinates
      _fe_u.coorMap(x, y, z, _fe_u.refFE.xi(lDof-1), _fe_u.refFE.eta(lDof-1), _fe_u.refFE.zeta(lDof-1));
	      
      // Loop on data vector components
      for (UInt icmp=0; icmp < nbComp; ++icmp)
	_u.vec()( icmp*_dim_u + _dof_u.localToGlobal(iElem,lDof) - 1) =   u0(0.0,x,y,z,icmp+1);      
    }
  }
}


#endif
