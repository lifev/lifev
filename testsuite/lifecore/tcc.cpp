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
#include<iostream>
class b;
class a;
class common
{
public:
typedef float Real;
typedef a A;
typedef b B; 
};

class a 
{
public:
  a(): a1(1), a2(10.){};
  a(int i, float f): a1(i), a2(f){}; 
  a(const a & s)
    {a1=s.g_a1();
    a2=s.g_a2();
      };
  int g_a1() const {return a1;}
  float g_a2() const {return a2;}
  a & operator =(const a & gg){a1=0;a2=gg.g_a2();}
private:
  int a1;
  float a2;
};
class b
{
public:
  b(): ib(90){};
  int ib;
};
template <class T>
class c: a , public T::B
{
public:
  friend 
  ostream & operator << (ostream & s, c<T> & x);
  typedef typename T::B b;
  c(){};
  c(int i , float f):a(i,f), b() {};
  c(const c<T> & z): a(z), b(z){};
  c<T> & operator =(const c<T> & gg) {a::operator=(gg);};  
};

template<class T >
ostream & operator << (ostream & s, c<T> & x);

template<>
ostream & operator <<  (ostream & s, c<common> & x)
{
  return s << x.ib<<", "<< x.g_a1() << ", "<<x.g_a2() << endl;
}

int
main()
{
  c<common> pippo(1,89.0);
  cout <<pippo;
  c<common> pluto;
  cout <<pluto;
  c<common> paperino(pippo);
  cout<<paperino;
  paperino=pluto;
  cout<<paperino;

}

