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
#include "subDomainNeighbors.hpp"

//! useful function to sort a list and remove multiple numbers.
void RemoveMultiple(const vector<ID> & list0, vector<ID> & listf){

  Index_t counter;
  vector<ID> tmplist = list0;

  //! Sort the list
  sort(tmplist.begin(), tmplist.end());
  
  //! initialize the new list
  listf.push_back( tmplist[0] );
  counter = 0;

  //! We remove the multiple occurences :
  for (Index_t i = 1 ;  i < tmplist.size() ; i++){
    if ( tmplist[ i ] != listf[ counter ] ){
      //! Add to the list the new value
      listf.push_back( tmplist[i] );
      counter ++ ;
    }
  }
}

//! Constructor
SubDomainNeighbors::SubDomainNeighbors(ID SDomID):
  _SDomID(SDomID),
  _nbInterf(0),
  _nbNeigh(0),
  _finalized(false)
{}

//! Constructor taking the connectivity table in a file as input (to do)
SubDomainNeighbors::SubDomainNeighbors(ID SDomID, string fname):
  _SDomID(SDomID),
  _nbInterf(0),
  _nbNeigh(0),
  _finalized(false)
{}

//! Destructor (useless?)
SubDomainNeighbors::~SubDomainNeighbors()
{}

//! How many neighbors stored?
Index_t SubDomainNeighbors::sizeNeigh() const{
  return _neighList.size(); 
}
  
//! return the reference of the interface i in the neighbors' list. (Beware: i starts from 0). 
Int SubDomainNeighbors::NeighInterfaceRef( const Index_t & i ) const{
  return _neighList[ i ].InterfaceRef; 
}


//! Output
ostream& SubDomainNeighbors::showMe(bool verbose, ostream & out) const {
  out << "********************************"<<endl;
  out << "SubDomain number: " << _SDomID << endl;  
  out << "Number of Neighbors: " << _nbNeigh << endl;
  out << "Size of _neighList: " << _neighList.size() << endl;
  out << "List of Neighbors: \n";
   // starts from 0
  for (Index_t i = 0 ;  i < _nbNeigh ; i++){
    out << _neighList[i].NeighborID <<" ";
    out << _neighList[i].InterfaceRef <<" ; ";
  }
  out << endl;
  out << "Number of Interfaces: " << _nbInterf << endl;
  out << "Size of _interfList: " << _interfList.size() << endl;
  out << "List of Interfaces: \n";
   // starts from 0
  for (Index_t i = 0 ;  i < _nbInterf ; i++)
    out << _interfList[i] <<" ";
  out << endl;

  out << "********************************" << endl;
  return out;
}
