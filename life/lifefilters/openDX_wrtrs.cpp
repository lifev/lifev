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
/*! --------------------------------------------------------------------------*
/                                                                            /
/      ...                                                                   /
/                                                                            /
/                                                                            /
/ OPENDX WRITERS UTILITIES                                                   /
/                                                                            /
/ #Version 0.1 Experimental:   30/08/2001 Alessandro Veneziani               /
/                                                                            /
/  I am sorry for my C++ at a very beginner level !!!!                       /
/                                                                            /
/ #Purpose: Container for subroutines for OPENDX                             /
/                                                                            /
/                                                                            /
/---------------------------------------------------------------------------*/
#include "lifeV.hpp"
#include <vector>
#include "openDX_wrtrs.hpp"

namespace LifeV
{

void wr_opendx_vector(std::string fname, std::string name, std::vector<Real> U,
                      std::vector<Real> V, std::vector<Real> W)
{
 unsigned int i;
 char quotes='"';
 std::ofstream ofile(fname.c_str(),std::ios::app);
 std::string Ext = ".data";
 std::string nameExt = name + Ext;

 ASSERT(ofile,"Error: Output file cannot be open");


  ofile << "object " << quotes << nameExt << quotes << std::endl;
  ofile << "class array ";
  ofile << "type double ";
  ofile << "rank " << VRANK << " shape " << NDIM;
  ofile << " items " <<  U.size() << " ";
  ofile << "data follows" << std::endl;

  for (i=0;i<U.size();++i){
   if (fabs(U[i])<EPS_DX) U[i]=0. ; // needed due to dx bugs
   if (fabs(V[i])<EPS_DX) V[i]=0. ; // needed due to dx bugs
   if (fabs(W[i])<EPS_DX) W[i]=0. ; // needed due to dx bugs
   ofile << U[i] << " " << V[i] << " " << W[i] << std::endl;
  }

  ofile << "attribute " << quotes << "dep" << quotes << " string " << quotes
        << "positions" << quotes << std::endl;
  ofile << std::endl;
  ofile << "object " << quotes << name << quotes << " class field"
        << std::endl;
  ofile << "component " << quotes << "positions" << quotes << " value "
        << quotes << "pos" << quotes << std::endl;
  ofile << "component " << quotes << "connections" << quotes << " value "
        << quotes << "con" << quotes << std::endl;
  ofile << "component " << quotes << "data" << quotes << " value " << quotes
        << nameExt << quotes << std::endl;
  ofile << std::endl;
  ofile << "end";
  ofile << std::endl;
}
}
