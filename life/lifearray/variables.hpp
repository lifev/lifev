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
#ifndef _VARIABLES_H
#define _VARIABLES_H
#include<vector>
#include<string>
#include<typeinfo>
#include "tab.hpp"

namespace LifeV
{
template<typename UnkType,int size,int nblock>
class Unknown {

 public:
  /*  Unknown(std::string name):_name(name),_size(size),_nb(nblock)
   {if (size==1) {unk=0}
    else (unk(size,nblock));
    return;};

  Unknown():_name(" "),_size(size)
   {if (size==1) {unk=0;}
    else unk(size,nblock);
    return;};
  */

  int size() {return _size;}; //physical dimension of the unknown
  int ndof() {return _ndof;}; //degrees of freedom for the unknown
  int nb() {return _nb;};
  std::string name() {return _name;};

  std::vector<UnkType> unk;
  void set_vect(int ndof){unk(ndof);};


  //bc_selection & c.;


 private:
 std::string _name;
 int _size;
 int _nb;
 int _ndof;
}



typedef Unknown<double,1,1> ScaUnk; // Generic Scalar unknown

typedef Unknown<ElemVec,nDimension,1> PhyUnk; // Physical Vector unknown (e.g. velocity in NS problem)

}
#endif
