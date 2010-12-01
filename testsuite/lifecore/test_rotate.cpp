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
#include<iostream>
#include<vector>
#include<algorithm>

int main()
{
    std::vector<int> a;
    a.resize(10,int(1));
    for (std::vector<int>::iterator i=a.begin()+1; i != a.end(); ++i)
    {
        *i += *(i-1);
    }
    for (std::vector<int>::iterator i=a.begin(); i != a.end(); ++i)
    {
        std::cout<< *i <<", ";
    }
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::rotate(a.begin(),a.end()-1,a.end());
    for (std::vector<int>::iterator i=a.begin(); i != a.end(); ++i)
    {
        std::cout<< *i <<", ";
    }
    std::cout<<std::endl;
}
