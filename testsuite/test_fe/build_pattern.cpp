/* -*- mode: c++ -*-
   This program is part of the LifeV library 
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
#include <iostream>
#include <utility>

/*
  a stupid code that generates automatically  the pattern of a
  "standard" finite element (just cut and past the result in the file
  feDef.h, and don't forget to delete the last "," ).
                                               Jean-Fred Gerbeau 11/99.
*/
int main()
{
  int nbNode  = 10;
  int nbPattern = nbNode*nbNode;
  int nbDiag  = nbNode;
  int nbUpper = nbNode*(nbNode-1)/2;
  int ip,i,j;
  std::pair<int,int> _pattern[nbPattern];
  // first : diagonal terms
  for(ip=1;ip<=nbDiag;ip++)
    _pattern[ip-1].first = _pattern[ip-1].second = ip;
  // second : upper terms and third : lower terms
  ip = nbDiag+1;
  for(i=1;i<nbNode;i++){
    for(j=i+1;j<=nbNode;j++){
      _pattern[ip-1].first   = _pattern[ip-1+nbUpper].second = i;
      _pattern[ip-1].second  = _pattern[ip-1+nbUpper].first  = j;
      ip++;
    }
  }

  std::cout << "// First : rows " << std::endl;
  std::cout << "static int _pf[" << nbPattern << "] = {";
  for(ip=0;ip<nbDiag;ip++) std::cout << _pattern[ip].first << "," ;
  std::cout << "// diag  entries" << std::endl;
  for(ip=nbDiag;ip<nbDiag+nbUpper;ip++) std::cout << _pattern[ip].first << "," ;
  std::cout << "// upper  entries" << std::endl;
  for(ip=nbDiag+nbUpper;ip<nbPattern;ip++)std::cout << _pattern[ip].first << "," ;
  std::cout << "}; // lower  entries" << std::endl << std::endl;

  std::cout << "// Second : column " << std::endl;
  std::cout << "static int _ps[" << nbPattern << "] = {";
  for(ip=0;ip<nbDiag;ip++) std::cout << _pattern[ip].second << "," ;
  std::cout << "// diag  entries" << std::endl;
  for(ip=nbDiag;ip<nbDiag+nbUpper;ip++) std::cout << _pattern[ip].second << "," ;
  std::cout << "// upper  entries" << std::endl;
  for(ip=nbDiag+nbUpper;ip<nbPattern;ip++)std::cout << _pattern[ip].second << "," ;
  std::cout << "}; // lower  entries" << std::endl;

}

