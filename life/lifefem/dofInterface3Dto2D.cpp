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
#include "dofInterface3Dto2D.hpp"

namespace LifeV
{
//! Constructor for interfacing Dof of the same type (RefFE)
/*!
  \param refFe the part of the reference FE that contains the dof patterns (nbDofPerEdge...)
  \param dof1 the Dof object of the mesh in which we want to make the computations
*/
DofInterface3Dto2D::DofInterface3Dto2D( const LocalDofPattern& refFE, const Dof& dof1):
  _interfRef(0),_refFE1(&refFE),_dof1(&dof1)
{
  _finalized = false;
}

//! Returns the reference of the interface
EntityFlag DofInterface3Dto2D::InterfaceRef() const{
  return _interfRef;
}

//! Returns the identity of the i-th elements in the (finalised) face list
//! (counting from 0 ' a la C')
ID DofInterface3Dto2D::operator[](const UInt& i) const {
  ASSERT_PRE(_finalized, "The face List should be finalised before being accessed");
  ASSERT_BD(i >= 0 && i < _faceList.size());
  return _faceList[i].first;  // _faceList must be a vector!
}


//! Assignment operator (we have a vector of DofInterface3Dto2D)
DofInterface3Dto2D & DofInterface3Dto2D::operator=(const DofInterface3Dto2D& dofi){
  _interfRef = dofi._interfRef;
  _refFE1 = dofi._refFE1;
  _dof1   = dofi._dof1;
  _faceList          = dofi._faceList;
  _vertexPerFaceList = dofi._vertexPerFaceList; // (empty)
  _vertexList        = dofi._vertexList;
  // _edgePerFaceList   = dofi._edgePerFaceList; // (empty)
  // _edgeList          = dofi._edgeList;
  _locDof            = dofi._locDof;            // (empty)
  _locDofMap         = dofi._locDofMap;
  _finalized         = dofi._finalized;

  return *this;
}

//! true if the lists have been updated.
bool DofInterface3Dto2D::finalized() const {
  return _finalized;
}

//! removes all unuseful list (all except _faceList). use it properly!
void DofInterface3Dto2D::ClearLists(){
  _vertexPerFaceList.clear();
  _vertexList.clear();
  //  _edgePerFaceList.clear();
  //  _edgeList.clear();
  _locDof.clear();
}


//! Transforms the 3d index of a vertex into its 2d (interface) index.
//! This is a simple algorithm... Find out something better some day...?
ID  DofInterface3Dto2D::_Vtx3Dto2D( const ID& idpoint3D ) const {
  ASSERT_PRE(_finalized, "The list of vertices must be finalized before accessing to the interface vertices." );
  for (std::list< std::pair<ID,ID> >::const_iterator it = _vertexList.begin(); it!=_vertexList.end(); ++it) {
    if ( it->first == idpoint3D )
      return  it->second;
  }
  ERROR_MSG("There is no such 3D index of vertex in the _vertexList.");
}


//! Output
std::ostream& DofInterface3Dto2D::showMe2D(bool verbose, std::ostream& out) const  {
  out << "------------------------------"<< std::endl;
  out << "myDofInterface reference: " << _interfRef  << std::endl;
  out << "Number of face connections (_faceList): " << _faceList.size() << std::endl;
  if ( verbose ){
    unsigned int count(0),lines(10);
    out << "\tList of connections between Faces: (global, local)";
    for (std::vector< std::pair<ID,ID> >::const_iterator i=_faceList.begin(); i != _faceList.end(); ++i) {
      if (count++ % lines ==0){
	out << std::endl;
      }
      out << "(" << i->first << "," <<  i->second << ")\t";
    }
    out << std::endl;
  }
  out << "Number of connections between Vertices (_vertexList): " <<  _vertexList.size() << std::endl;
  if ( verbose ){
    unsigned int count(0),lines(10);
    out << "\tList of connections between Vertices: (global, local)";
    for (std::list< std::pair<ID,ID> >::const_iterator it = _vertexList.begin(); it!=_vertexList.end(); ++it) {
      if (count++ % lines ==0){
	out << std::endl;
      }
      out << "(" << it->first << "," << it->second << ")\t";
    }
    out << std::endl;
  }
  //! print _locDofMap
  showMe(verbose, out);

  out << "------------------------------" << std::endl;
  return out;
}



//! useful function to sort a list and remove multiple numbers.
void RemoveMultiple(const std::list<ID> & list0, std::list< std::pair<ID,ID> > & listf){

  ID counter = 1;
  std::list<ID> tmplist(list0);

  //! Sort the list
  tmplist.sort();

  //! initialize the new list
  std::pair <ID,ID>  p0( tmplist.front() , counter );
  listf.push_back( p0 );

  //! We remove the multiple occurences :
  for (std::list<ID>::iterator it = tmplist.begin() ;  it !=tmplist.end() ; ++ it ){
    if ( (*it) != listf.back().first ){
      counter ++ ;
      //! Add to the list the new value
      std::pair <ID,ID>  p( (*it) , counter );
      listf.push_back( p );
    }
  }
}
}
