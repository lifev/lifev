/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003 LifeV Team
  
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
  \file dof.h
  \brief Degrees of freedom, the class that provides the localtoglobal table
  \version 1.0
  \author M.A. Fernandez & Luca Formaggia
  \date 07/2002

  The numbering of the degrees of freedom follow the order
  Vertices - Edges - Faces - Volumes
  If more than one degree of freedon is associated to a geometrical entity, the
  global numbering associated to two "consecutive" degree of freedom  on the given
  entity is consecutive, and the order is according to the entity orientation.


  \author Modified by Vincent MARTIN & Mohamed BELHADJ
  \date 19/07/2002
   
  We added the HdivFE element here. 
*/

#ifndef _DOF_HH
#define _DOF_HH

#include "lifeV.hpp"
#include "SimpleVect.hpp"
#include "localDofPattern.hpp"
#include <algorithm>

/*! Local-to-global table

 This class provides the localtoglobal table that relates the local DOF of
 a finite element to its global numbering. It needs a LocalDofPattern in
 order to obtain all the necessary information about the local pattern. In
 fact it stores a copy of it so to make the local pattern available, if
 needed.

 It is useless until is has not been set up on a specific RegionMesh. This is accomplished either by
 passing the mesh to the constructor, or calling the method Dof::update().

 \note The methods bulds the table for ALL degrees of freedom, i.e. it does not handle any essential
 boundary condition.
*/
class Dof {
public:
  //! Type for the localToGlobal table.
  typedef SimpleArray<UInt> Container;
  
  //! The pattern of the local degrees of freedom.
  /*! It is exposed so that is is possible to interrogate it directly.*/
  const LocalDofPattern& fe; // be careful : use fe.nbLocalDof (fe.nbDof does not exist !)
  
  /*! The minimal constructor
    \param _fe is the LocalDofPattern on which the ref FE is built
    \param Offset: the smallest Dof numbering. It might be used if we want the
    degrees of freedom numbering start from a specific value.
  */
  Dof(const LocalDofPattern& _fe, UInt offSet=1);

  Dof(const Dof & dof2);
  
  //! Constructor accepting a mesh as parameter
  /*!
    \param mesh a RegionMesh3D
    \param _fe is the LocalDofPattern on which the ref FE is built
    \param Offset: the smalest Dof numbering. It might be used if we want the
    degrees of freedom numbering start from a specific value.
  */
  template<typename Mesh> 
    Dof(Mesh& mesh, const LocalDofPattern& _fe, UInt offSet=1);

  //! Build the localToGlobal table
  /*!
    \param mesh A RegionMesh3D
    Updates the LocaltoGlobal array
  */
  template <typename Mesh> 
    void update(Mesh &);
  
  //! The total number of Dof
  inline UInt numTotalDof() const {return _totalDof;}

  //! The number of local Dof (nodes) in the finite element
  inline UInt numLocalDof() const {return fe.nbLocalDof;}

  //! Return the specified entries of the localToGlobal table
  /*!
    Returns the global numbering of a DOF, given an element and the local numbering
    \param ELId the element ID
    \param localNode the local DOF numbering (starting from 1)
    \return The numbering of the DOF
  */
  inline ID localToGlobal(const ID ElId, const ID localNode) const {
    return _ltg(localNode,ElId);
  }

  //! Number of elements in mesh 
  UInt numElements() const {return _nEl;}
  
  //! Number of local vertices (in a elment)
  UInt numLocalVertices() const {return nlv;}
  
  //! Number of local edges (in a elment)
  UInt numLocalEdges() const {return nle;}
  
  //! Number of local faces (in a elment)
  UInt numLocalFaces() const {return nlf;}
  
  //! Ouput
  void showMe(std::ostream  & out=std::cout, bool verbose=false) const;

private:
  UInt _offset;
  UInt _totalDof;
  UInt _nEl;  
  UInt nlv;
  UInt nle;
  UInt nlf;
  Container _ltg;
  UInt _ncount[5];
};
 


/********************************************************************
                      IMPLEMENTATIONS
********************************************************************/

//! Constructor that builds the localToglobal table
template <typename Mesh> 
Dof::Dof(Mesh& mesh, const LocalDofPattern& _fe, UInt off):fe(_fe),_offset(off),_totalDof(0),
					       _nEl(0),nlv(0),nle(0),nlf(0),_ltg()
{ for (UInt i=0; i<5; ++i)_ncount[i]=0;
  update(mesh);
}


//! Build the localToGlobal table
template<typename Mesh> 
void Dof::update(Mesh& M){

  typedef  typename Mesh::VolumeShape GeoShape;

  // Some useful local variables, to save some typing
  UInt nldpe=fe.nbDofPerEdge;
  UInt nldpv=fe.nbDofPerVertex;
  UInt nldpf=fe.nbDofPerFace;
  UInt nldpV=fe.nbDofPerVolume;

  nlv=GeoShape::numVertices;
  nle=GeoShape::numEdges;
  nlf=GeoShape::numFaces;

  _nEl=M.numElements();

  UInt nV=M.numVolumes();
  UInt ne=M.numEdges();
  UInt nv=M.numVertices();
  UInt nf=M.numFaces();

  UInt i,l,ie;
  
  UInt nldof=nldpV+nldpe*nle+nldpv*nlv+nldpf*nlf;

  ASSERT_PRE( nldof == UInt(fe.nbLocalDof), "Something wrong in FE specification") ;
  
  _totalDof=nV*nldpV+ne*nldpe+nv*nldpv+nf*nldpf;
  
  _ltg.reshape(nldof,nV);
  
  // Make sure the mesh has everything needed
  bool update_edges(nldpe !=0 && ! M.hasLocalEdges());
  bool update_faces(nldpf !=0 && ! M.hasLocalFaces());
  
  if (update_edges) M.updateElementEdges();
  if (update_faces) M.updateElementFaces();
  //  ASSERT_PRE( !(nldpe !=0 && M.hasLocalEdges()) , "Element edges stuff have not been updated") ;
  //  ASSERT_PRE( !(nldpf !=0 && M.hasLocalFaces()) , "Element faces stuff have not been updated") ;
  //ASSERT_PRE( (nldpe == 0 || M.hasLocalEdges()) , "Element edges stuff have not been updated") ;
  //ASSERT_PRE( (nldpf == 0 || M.hasLocalFaces()) , "Element faces stuff have not been updated") ;
  

  unsigned int gcount(_offset);
  unsigned int lcount;
  unsigned int lc;
  
  // Vertex Based Dof
  _ncount[0]=gcount;
  if (nldpv >0 )
    for (ie=1; ie<= _nEl; ++ie){
      lc=0;
      for (i=1; i<=nlv; ++i)
	for (l=0; l<nldpv; ++l)
	  _ltg(++lc,ie)=gcount+(M.volume(ie).point(i).id()-1)*nldpv + l;
    }
  // Edge Based Dof
  gcount+=nldpv*nv;
  lcount=nldpv*nlv;
  _ncount[1]=gcount;
  if (nldpe >0 )
    for (ie=1; ie<= _nEl; ++ie){
      lc=lcount;
      for (i=1; i<=nle; ++i)
	for (l=0; l<nldpe; ++l)
	  _ltg(++lc,ie)=gcount +(M.localEdgeId(ie,i)-1)*nldpe + l;
    }
 
  // Face  Based Dof
  
  gcount+=ne*nldpe;
  lcount+=nldpe*nle;
  _ncount[2]=gcount;
  if (nldpf >0 )
    for (ie=1; ie<= _nEl; ++ie){
      lc=lcount;
      for (i=1; i<=nlf; ++i)
	for (l=0; l<nldpf; ++l)
	  _ltg(++lc,ie)=gcount +(M.localFaceId(ie,i)-1)*nldpf + l;
    }

  // Volume  Based Dof
  gcount+=nf*nldpf;
  lcount+=nldpf*nlf;
  _ncount[3]=gcount;
  if (nldpV >0 )
    for (ie=1; ie<= _nEl; ++ie) {
      lc=lcount;
      for (l=0; l<nldpV; ++l)
	_ltg(++lc,ie)=gcount+(ie-1)*nldpV +l;
    }
  gcount+=nV*nldpV;
  _ncount[4]=gcount;
  ASSERT_POS(gcount -_offset == _totalDof , "Something wrong in Dof Setup "<<gcount<<" " <<_offset <<" " <<_totalDof) ;
  if (update_edges) M.cleanElementEdges();
  if (update_faces) M.cleanElementFaces();
}

#endif
