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
#include "geo0D.hpp"

/*--------------------------------------------------------------
                                 Geo0D
  ---------------------------------------------------------------*/
Geo0D::Geo0D():
MeshEntityWithBoundary(0)
{}

Geo0D::Geo0D(ID id,bool boundary):
MeshEntityWithBoundary(id,boundary)
{}

Geo0D::Geo0D(ID id, Real x, Real y, Real z, bool boundary):
MeshEntityWithBoundary(id,boundary)
{
_coor[0]=x;
_coor[1]=y;
_coor[2]=z;
}

Geo0D::Geo0D(Geo0D const & G):
MeshEntityWithBoundary(G._id,G._boundary)
{
  for (UInt i =0; i< nDimensions; ++i)
    _coor[i]=G._coor[i];
}

Geo0D &
Geo0D::operator=(Geo0D const & G)
   // Assignement operator
{
  _id=G._id;
  _boundary=G._boundary;
  for (UInt i =0; i< nDimensions; ++i)
    {
    _coor[i]=G._coor[i];
    }
  return *this;
}


ostream & Geo0D::showMe(bool verbose, ostream & out) const
{
  out.setf(ios::scientific,ios::floatfield);
  out << " Geo0D object " <<endl;
  if (verbose) 
    {
      out << " Coordinates:" <<endl;
      Real const * c=coor();
      for ( unsigned i=0; i<nDimensions; i++ )
	{
	  out << c[i] <<",  ";
	}
      out << endl<<endl;
    }
  out << "ID= "<< id()<<"  ";
  out << "----- END OF Geo0D data ---"<<endl<<endl;
  return out;
};

