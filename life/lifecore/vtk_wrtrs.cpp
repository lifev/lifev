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
#include "vtk_wrtrs.hpp"


void wr_vtk_ascii_scalar(string fname, string name, Real* U, int Usize,
			 string look_up_table)
{
 int i,j;
 ofstream ofile(fname.c_str(),ios::app);

 ASSERT(ofile,"Error: Output file cannot be open");
 
 ofile << "SCALARS " << name << " " << "float " << "1" << endl;
 ofile << "LOOKUP_TABLE " << look_up_table << endl;

 for (i=0;i<Usize;++i)
  ofile << U[i] << endl;

 if (look_up_table!="default")
   {
    ofile << "LOOKUP_TABLE " << look_up_table << " " << Usize << endl; 
    ifstream lutfile(look_up_table.c_str());
    ASSERT(lutfile,"Error: LookUp Table  file cannot be open"); 
    // assert
    Real r,g,b,a;
    for (j=0;j<Usize;++j){
     lutfile >> r >> g >> b >> a;
     ofile << r << " " << g << " " << b << " " << a << endl;   
    }
   } 
}
void wr_vtk_ascii_vector(string fname, string name, Real* U, int Usize)
{
 int i,j;
 ofstream ofile(fname.c_str(),ios::app);
 int nbcomp = 3;

 ASSERT(ofile,"Error: Output file cannot be open"); 


 ofile << "VECTORS " << name << " " << "float " << endl;
 for (i=0;i< Usize/nbcomp; ++i){
     for (j=0; j<nbcomp; j++)
       ofile << U[i+j*Usize/nbcomp] << " ";
     ofile << endl;
   }
}

//----------------------------------------------------------------------
// obsolete ?
void wr_vtk_ascii_scalar(string fname, string name, vector<Real> U, string look_up_table)
{
 unsigned int i,j;
 ofstream ofile(fname.c_str(),ios::app);

 ASSERT(ofile,"Error: Output file cannot be open"); 


 ofile << "POINT_DATA " << U.size() << endl;
 ofile << "SCALARS " << name << " " << "float " << "1" << endl;
 ofile << "LOOKUP_TABLE " << look_up_table << endl;

 for (i=0;i<U.size();++i)
  ofile << U[i] << endl;

 if (look_up_table!="default")
   {
    ofile << "LOOKUP_TABLE " << look_up_table << " " << U.size() << endl; 
    ifstream lutfile(look_up_table.c_str());
    ASSERT(lutfile,"Error: LookUp Table  file cannot be open"); 
    // assert
    Real r,g,b,a;
    for (j=0;j<U.size();++j){
     lutfile >> r >> g >> b >> a;
     ofile << r << " " << g << " " << b << " " << a << endl;   
    }
   } 
}

