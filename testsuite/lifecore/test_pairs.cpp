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
#include <utility>
#include <iostream>
int
main()
{
    std::pair<int,int> a,b,c;
    int i1=1;
    int i2=3;
    int i3=5;
    a=std::make_pair(i1,i2);
    b=std::make_pair(i2,i1);
    std::cout << " a less than b? "<< (a < b) << std::endl;
    std::pair<int, std::pair<int,int> > a1,a2;
    a1=std::make_pair(i1,std::make_pair(i2,i3));
    a2=std::make_pair(i1,std::make_pair(i3,i2));
    std::cout << " a1 less than a2? "<< (a1 < a2) << std::endl;
}

