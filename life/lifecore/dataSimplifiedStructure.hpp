/*!
  \file dataSimplifiedStructure.h
  \author M.A. Fernandez
  \date 01/2003 
  \version 1.0

  \brief File containing a class for handling data for reduced 
         structural models (algegraic law and independent ring model)

*/

#ifndef _DATASIMPLIFIEDSTRUCTURE_H_
#define _DATASIMPLIFIEDSTRUCTURE_H_
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "lifeV.hpp"

/*! 
  \class DataNavierStokes

  Base class which holds usual data for reduced 
  structural models (algegraic law and independent ring model)

*/
template<typename Mesh>
class DataSimplifiedStructure:
public DataMesh<Mesh> {
 public:

  //! Constructor
  DataSimplifiedStructure(const GetPot& dfile);
  
  //! Ouptut
  void showMe(ostream& c=cout) const;

 protected:
  //! Physics
  Real _rho; // densisty
  Real _h;  // thickness
  Real _E;  // Young modulus
  Real _nu; // Poisson coeficient
  Real _R0; // Radius
};



//
// IMPLEMENTATION
//


// Constructor
template<typename Mesh>
DataSimplifiedStructure<Mesh>::
DataSimplifiedStructure(const GetPot& dfile):
  DataMesh<Mesh>(dfile,"solid/discretization")
{
  // physics
  _rho     = dfile("solid/physics/density",1.);
  _h       = dfile("solid/physics/thickness",1.);
  _E       = dfile("solid/physics/young",1.);
  _nu      = dfile("solid/physics/poisson",0.5);
  _R0      = dfile("solid/physics/radius",1.);
}


// Output
template<typename Mesh>
void DataSimplifiedStructure<Mesh>::
showMe(ostream& c) const
{
  // physics
  c << "\n*** Values for data [solid/physics]\n\n";
  c << "density   = " << _rho << endl; 
  c << "thickness = " << _h   << endl;
  c << "young     = " << _E   << endl; 
  c << "poisson   = " << _nu  << endl; 
  c << "radius    = " << _R0  << endl; 
  c << "\n*** Values for data [solid/discretization]\n\n";
  DataMesh<Mesh>::showMe();
}

#endif
