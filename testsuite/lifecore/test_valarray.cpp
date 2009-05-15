/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano
  
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
#include <iostream>
#include <valarray>

template <typename T>
T dot(std::valarray<T> const a, std::valarray<T> const b)
{
  T c=0;
  unsigned int  sm;
  sm = a.size() < b.size() ? a.size(): b.size();

  for (unsigned int i=0; i< sm; ++i){
    c+=a[i]*b[i];
  }
  return c;
  
}

int
main()
{
  std::valarray<float> f1(100);
  std::valarray<float> f2(100);
  std::valarray<int> i1(100);
  f1=50.0;
  for (unsigned int i=0; i< f2.size(); i++){
    
    f2[i]=i+0.5;
  }
  for (unsigned int i=0; i< i1.size(); i++){
    
    i1[i]=i;
  }
  std::cout << f2[5]<< " " <<f1[49] <<std::endl;
  std::cout << i1[90]<<std::endl;
  std::cout << i1.sum()<<" " << i1.size()*(i1.size()-1)*0.5<<std::endl;

  i1+=10;
  std::cout << i1.sum()<<std::endl;

  f1[2]=40.5;
  f1=f1+f2;
  std::cout << f1[2]<<std::endl;
  std::cout << dot(f1,f2)<<std::endl;
  
  
}

