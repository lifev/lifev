/*!
  \file dataElasticStructure.h
  \author M.A. Fernandez
  \date 10/2003 
  \version 1.0

  \brief 


*/

#ifndef _DATAELASTICSTRUCTURE_H_
#define _DATAELASTICSTRUCTURE_H_
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "lifeV.hpp"
#include "dataMesh.hpp"
#include "dataTime.hpp"
/*! 
  \class DataElasticStructure

 

*/
template<typename Mesh>
class DataElasticStructure:
public DataMesh<Mesh>, 
public DataTime {
 public:

  //! Constructor
  DataElasticStructure(const GetPot& dfile);
  
  //! Ouptut
  void showMe(ostream& c=cout) const;

  //! End time
  Real endtime() const;

 protected:
  //! Physics
  Real _rho; // densisty
  Real _E;  // Young modulus
  Real _nu; // Poisson coeficient
  Real _lambda, _mu; // Lame coefficients
  Real _endtime; // end time

  //! Miscellaneous  
  Real _factor; // amplification factor for deformed mesh
  UInt _verbose; // temporal output verbose

};



//
// IMPLEMENTATION
//


// Constructor
template<typename Mesh>
DataElasticStructure<Mesh>::
DataElasticStructure(const GetPot& dfile):
  DataMesh<Mesh>(dfile,"solid/discretization"),
  DataTime(dfile,"solid/discretization") {
  // physics 
  _rho       = dfile("solid/physics/density",1.);
  _E       = dfile("solid/physics/young"    ,1.);
  _nu      = dfile("solid/physics/poisson"   ,0.25);
  _endtime   =  dfile("solid/physics/endtime",1.);

  // miscellaneous
  _factor    = dfile("solid/miscellaneous/factor",1.0);
  _verbose   = dfile("solid/miscellaneous/verbose",1);
  
  // Lame coefficients
  _lambda = _E*_nu / ( (1.0+_nu)*(1.0-2.0*_nu) );
  _mu     = _E /( 2.0*(1.0+_nu) );

}


// Output
template<typename Mesh>
void DataElasticStructure<Mesh>::
showMe(ostream& c) const
{
  // physics
  c << "\n*** Values for data [solid/physics]\n\n";
  c << "density                          = " << _rho << endl; 
  c << "young                            = " << _E   << endl; 
  c << "poisson                          = " << _nu  << endl; 
  c << "lame constants (lambda, mu)      = " << _lambda << " " << _mu  << endl;  
  c << "endtime                          = " << _endtime << endl;

  c << "\n*** Values for data [solid/miscellaneous]\n\n";
  c << "deformation factor               = " << _factor << endl;
  c << "verbose                          = " << _verbose << endl;

  c << "\n*** Values for data [solid/discretization]\n\n";
  DataMesh<Mesh>::showMe();
  DataTime::showMe(c);
}


// The end time
template <typename Mesh>
Real DataElasticStructure<Mesh>::
endtime() const {
  return  _endtime;
}


#endif
