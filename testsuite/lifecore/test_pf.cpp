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

class pluto
{
public:
    static int p1(int i)
    { return i;}
    static int p2(int i)
    { return i*2;}
};

template<class T>
class prova
{
public:
    typedef int (*Pfun)(int );
    Pfun pippo(int i)
    {
        if (i ==1 )
        {
            return T::p1;
        }
        else
        {
            return T::p2;
        }
    }
};

int main ()
{
    prova<pluto> a;
    std::cout<<a.pippo(1)(1)<<std::endl;
    std::cout<<a.pippo(2)(1)<<std::endl;
}
