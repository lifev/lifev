/*!
  \file dataConvDiffReact.hpp
  \author M. Prosi
  \date 03/2004 
  \version 1.0

  \brief File containing a class for handling Convection-Diffusion-Reaktion processes data with GetPot

*/
#ifndef _DATACONVDIFFREACT_H_
#define _DATACONVDIFFREACT_H_
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "lifeV.hpp"
#include "dataMesh.hpp"
#include "dataTime.hpp"

/*! 
  \class DataConvDiffReact

  Base class which holds usual data for the Convection-Diffusion-Reaction equation solvers

*/
template <typename Mesh>
class DataConvDiffReact:public DataMesh<Mesh>,public DataTime {
 public:

  //! Constructor
  DataConvDiffReact(const GetPot& dfile);
  
  //! Ouptut
  void showMe(ostream& c=cout);
  //! Diffusivity
  Real diffusivity() const;
  //! Reaction coefficient
  Real react() const;

 protected:
  //! Physics
  Real _diffusivity; // Diffusivity
  Real _react; // Reaction coefficient
};


//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh>
DataConvDiffReact<Mesh>::
DataConvDiffReact(const GetPot& dfile):
  DataMesh<Mesh>(dfile,"masstransport/discretization"),
  DataTime(dfile,"fluid/discretization") {
  
  // physics
  _diffusivity = dfile("masstransport/physics/diffusivity",1.);
  _react = dfile("masstransport/physics/react",1.);
}

// Output
template <typename Mesh>
void DataConvDiffReact<Mesh>::
showMe(ostream& c)
{
  // physics
  c << "\n*** Values for data [masstransport/physics]\n\n";
  c << "diffusivity   = " << _diffusivity << endl; 
  c << "reaction coefficient  = " << _react << endl;

  c << "\n*** Values for data [masstransport/discretization]\n\n";
  DataMesh<Mesh>::showMe(c);
  DataTime::showMe(c);

}
////////////////////
// The diffusivity
template <typename Mesh>
Real DataConvDiffReact<Mesh>::
diffusivity() const {
  return  _diffusivity;
}

// The reaction coefficient
template <typename Mesh>
Real DataConvDiffReact<Mesh>::
react() const {
  return  _react;
}

#endif
