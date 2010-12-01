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
#include <vector>
#include <string>
class Stupid
{
public:
    Stupid(Stupid const & c)
    {
        message=c.message;
        v=c.v;

        std::cout << "copy costructor "<< message <<std::endl;
    }

    Stupid()
    {
        message="da standard constructor";
        std::cout << "building "<< message <<std::endl;
        v.resize(10);
        v[9]=90;
    }

    Stupid(std::string mess)
    {
        message=mess;
        std::cout << "building "<< message <<std::endl;
        v.resize(10);
        v[9]=90;
    }
    ~Stupid()
    {
        std::cout << "destroying "<< message <<std::endl;
    }
    void
    showme()
    {
        std::cout << v[9] <<" " << std::endl;
    }

private:
    std::string message;
    std::vector<int> v;

};

Stupid
crea1()
{
    Stupid a("funzione passa valore");
    return a;
}

Stupid
&
crea2()
{
    Stupid * a= new Stupid("funzione passa ref");
    return *a;

}



void ruotine()
{
    std::cout << "assegno  valore a valore"<<std::endl;
    Stupid a=crea1();
    a.showme();

    std::cout << "assegno valore a ref"<<std::endl;
    Stupid b =crea1();
    b.showme();
    std::cout << "assegno ref a ref"<<std::endl;
    Stupid & c =crea2();
    c.showme();
    std::cout << "assegno a valore a ref"<<std::endl;
    Stupid d =crea2();
    c.showme();
    delete &c;;

}

int main()
{
    std::cout<<" Calling"<<std::endl;

    ruotine();
    std::cout<<"Returned"<<std::endl;

}
