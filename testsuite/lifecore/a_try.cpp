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
#include <typeinfo>
using namespace std;

template<class TTTT> ostream & typeName(TTTT const &, ostream & out);

template<class T>
ostream & typeName(T const & a, ostream & out){
  out<<typeid(a).name();
  return out;
};


class base
{
public:
  typedef unsigned int UInt;
  typedef float Real;
  float ggg;
  
};


class derived : public base
{
public:
  
  derived():base(),c(10)
    {}
  UInt c;
  Real d;
  
};

template <typename Pippo>
class pluto : public Pippo
{
public:
  typedef Pippo Pluto;
  
  int ff;
};



int
main()
{
  derived a;
  typedef derived::UInt k;
  pluto<base> cc;
  typedef pluto<base>::UInt kl;
  kl yy=30;

  typedef pluto<base>::Pluto BB;
  BB i;
  
  i.ggg=90.9;
  
  
  std::cout <<a.c<<" " <<yy<<" "<<i.ggg<<endl;
}
