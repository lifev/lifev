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
#include "dofInterfaceBase.hpp"


//! Default Constructor
DofInterfaceBase::DofInterfaceBase(){}


void DofInterfaceBase::ReadVectorDataAndDofMap(const string filename, Vector& dataVec){

  UInt vecsize, idof;
  Real val;

  ifstream ifile(filename.c_str());

  ASSERT(ifile,"Error: Input Dof External BC file cannot be opened.");
 
  ifile >> vecsize;
  ASSERT(  dataVec.size() == vecsize, 
          "The vector in input must have the same dimension as the interface vector.");

  for(UInt i = 0; i< vecsize; i++){   
    ifile >> idof >> val;
    dataVec[i] = val;
    _locDofMap[idof] = i+1;
  }
}


//! This method returns the corresponding dof number of the mesh2 at the interface 
//! for a specific dof number at the interface in mesh1
/*!
  \param i a dof number in mesh1 
*/
ID DofInterfaceBase::getInterfaceDof(const ID& i) const {
  map<ID,ID>::const_iterator it  = _locDofMap.find(i);
  if (it == _locDofMap.end()) 
    ERROR_MSG("Dof number not found");
  return it->second;
} 


//! This method returns the number of dof that live on the interface
ID DofInterfaceBase::nbInterfaceDof() const{
  return _locDofMap.size();
}

std::ostream& DofInterfaceBase::showMe(bool verbose, std::ostream& out) const  {
  out << "------------------------------"<< std::endl;
  out << "\tNumber of Dof connections (_locDofMap):" << _locDofMap.size() << std::endl;
  if ( verbose ){
    UInt count(0),lines(10);
    out << "List of connections between Dof: (global, local)"; 
    for (map<ID,ID>::const_iterator it = _locDofMap.begin(); it!=_locDofMap.end(); ++it) {
      if (count++ % lines ==0){
	out << std::endl;
      }
      out << "(" << it->first << "," << it->second << ")\t";
    }
    out << std::endl;  
  }
  out << "------------------------------" << std::endl;
  return out;
}
