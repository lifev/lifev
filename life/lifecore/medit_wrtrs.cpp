/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003 LifeV Team
  
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
#include "medit_wrtrs.hpp"
#include "currentFE.hpp"

void wr_medit_ascii_scalar(string fname, Real* U, int Usize,int type)
{
  ofstream ofile(fname.c_str());
  
  ASSERT(ofile,"Error: Output file cannot be opened.");
 
  ofile << nDimensions << " 1 " << Usize << " " << type << endl; 

  for(int i = 0; i< Usize; i++){   
    ofile << U[i] << endl;
  }
}

void wr_medit_ascii_vector(string fname, Real* U, int Usize,int type)
{
  ofstream ofile(fname.c_str());

  ASSERT(ofile,"Error: Output file cannot be opened.");
 
  ofile << nDimensions << " " << nDimensions << " "
	<< Usize/nDimensions << " "  << type << endl; 
  
  for (int i=0;i<(int)Usize/(int)nDimensions; i++){
    for (int j=0; j<(int)nDimensions; j++)
      ofile << U[i+j*Usize/nDimensions] << " ";
    ofile << endl;
  }
}


void rd_medit_ascii_scalar(string fname, Real* U, const UInt& Usize, UInt& type)
{

  UInt theDim,theSize,nCol;
  ifstream ifile(fname.c_str());

  ASSERT(ifile,"Error: Input file cannot be opened.");
 
  ifile >> theDim >> nCol >> theSize >> type; 

  ASSERT_PRE( theDim==3 && nCol==1, 
	      "You cannot use rd_medit_ascii_scalar for non medit files or non scalar data.");
  ASSERT( Usize == theSize, 
	     "The vector in input must have the same dimension as the interface vector.");

  for(UInt i = 0; i< Usize; i++){   
    ifile >> U[i];
  }
}

void rd_medit_ascii_vector(string fname, Real* U, const UInt& Usize, UInt& type)
{
  UInt theDim,theSize,nCol;
  ifstream ifile(fname.c_str());

  ASSERT(ifile,"Error: Intput file cannot be opened.");

  ifile >> theDim >> nCol >> theSize >> type; 

  ASSERT_PRE( theDim==3 && nCol==nDimensions, 
	      "You cannot use rd_medit_ascii_vector for non medit files or non vector data.");
  ASSERT( theSize == Usize/nDimensions, 
	     "The vector in input must have the same dimension as the interface vector.");

  for (UInt i = 0 ; i < Usize/nDimensions ; i++ ){
    for (UInt j = 0 ; j < nDimensions ; j++)
      ifile >> U[i+j*Usize/nDimensions] ;
  }
}


