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

