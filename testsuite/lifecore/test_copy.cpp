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
#include <iostream>
#include <vector>
#include <algorithm>
template< typename sequence >
void
mycopy(std::vector<int>::iterator start, std::vector<int>::iterator end, sequence & to)
{
    std::copy(start,end,to);
}
int main () {
  std::vector<int> a(10);
  int b[10];
  for (int i =1; i<10; ++i) a[i]=i;
  mycopy(a.begin(),a.end(),b);
  std::cout<< b[6] << "\n";
}

