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
/*!
  \file simplifiedStructure.h
  \author M.A. Fernandez
  \date 01/2003 
  \version 1.0

  \brief File containing a classes handling reduced 
         structural models (algegraic law and independent ring model)

*/

#ifndef HH_SIMPLIFIEDSTRUCTURE_HH_
#define HH_SIMPLIFIEDSTRUCTURE_HH_

#include "lifeV.hpp"
#include "markers.hpp" 
#include "vecUnknown.hpp"
#include "dof.hpp"
#include "dataSimplifiedStructure.hpp"

/*! 
  \class AlegbraicLaw

  Base class for reduced algegraic law structural model.

*/
template<typename Mesh>
class AlgebraicLaw:
public DataSimplifiedStructure<Mesh> {

 public:

  //! Constructor
  /*!
    \param dfile GetPot data file
    \param mesh marker for the structure 
  */
  AlgebraicLaw(const GetPot& dfile, const EntityFlag& marker);
  
  //! Compute dispalcement for a given pressure load
  void updateDisplacement(const Vector& p);

  //! Displacements at steps n, n-1 and n-2
  Vector& d_n()   {return _d_n.vec();  }
  Vector& d_n_1() {return _d_n_1.vec();}
  Vector& d_n_2() {return _d_n_2.vec();}

  //! Dof 
  Dof& dof() {return _dof;}

  //! Updates displacements at steps n, n-1 and n-2
  void timeAdvance( const Vector& d_n);

 private:
  
  Dof _dof; // dof
  UInt _dim_mesh; // number of total dof
  PhysVectUnknown<Vector> _d_n; // displacement at step n
  PhysVectUnknown<Vector> _d_n_1; // displacement at step n-1
  PhysVectUnknown<Vector> _d_n_2; // displacement at step n-2
  EntityFlag _marker; // mesh marker for the structure  
};


/*! 
  \class IndependentRing

  Base class which implements the independent ring structural model 
       using a mid-point rule time-stepping

*/
template <typename Mesh>
class IndependentRing:
public DataSimplifiedStructure<Mesh>,
public DataTime {
 public:

  //! Constructor 
  /*!
    \param dfile GetPot data file
    \param mesh marker for the structure 
  */
  IndependentRing(const GetPot& dfile, const EntityFlag& marker);
		
  //! Get initial condition
  void initialize();

  //! Updating right hand side
  void timeAdvance(const Vector& d_n);
  
  //! Updates displacement
  void iterate(const Vector& p);
  
  // Displacements at steps n, n-1 and n-2
  Vector& d_n()   {return _d_n.vec();  }
  Vector& d_n_1() {return _d_n_1.vec();}
  Vector& d_n_2() {return _d_n_2.vec();}

  //! Dof
  Dof& dof() {return _dof;}

  //! Relaxation of the displacement
  void relaxation(const Real& omega, Vector& d);

  //! Output
  void showMe(ostream& c=cout) const;

 private:

  EntityFlag _marker; // mesh marker for the structure

  Dof _dof; // dof

  Vector _eta_nplus1; // radial displacement 
  
  Vector _w_nplus1;  // radial velocity

  Vector _eta_nplus1old; // radial displacement (previous) 
  
  Vector _w_nplus1old; // radial displacement (prevoius)
  
  Vector _f; // Auxiliary vectors (right hand side displacement)
 
  Vector _g; // Auxiliary vectors (right hand side velocity)
  
  PhysVectUnknown<Vector> _d_n;   // 3D displacement at step n
  PhysVectUnknown<Vector> _d_n_1; // 3D displacement at step n-1
  PhysVectUnknown<Vector> _d_n_2; // 3D displacement at step n-2
  
  Real a, b, c, d; // Auxiliary constants
};


//
// IMPLEMENTATION
//

// Constructor
template<typename Mesh>
AlgebraicLaw<Mesh>::
AlgebraicLaw(const GetPot& dfile, const EntityFlag& marker):
  DataSimplifiedStructure<Mesh>(dfile),
  _dof(_mesh,_mesh.getRefFE()),
  _dim_mesh(_dof.numTotalDof()),
  _d_n(_dim_mesh),
  _d_n_1(_dim_mesh),
  _d_n_2(_dim_mesh) {

  _d_n.vec()  =0;
  _d_n_1.vec()=0.0; 
  _d_n_2.vec()=0.0;
  _marker = marker;
}

// Compute dispalcement for a given pressure load
template<typename Mesh> void AlgebraicLaw<Mesh>::
updateDisplacement(const Vector& p) {

  Real dR,R,x,y;
  
  for (UInt i=0; i < _dim_mesh; ++i) {
    
    if ( _mesh.point(i+1).marker() == _marker) {
      x = _mesh.point(i+1).x();    
      y = _mesh.point(i+1).y();    

      R  = sqrt(x*x+y*y);
      dR = (1.0-_nu*_nu)*_R0*_R0*p[i]/(_h*_E);
      
      _d_n.vec()[i]               = x * dR / R;
      _d_n.vec()[i + _dim_mesh]   = y * dR / R;
    }
  }  
}

 // Updates displacements at steps n, n-1 and n-2
template<typename Mesh> void AlgebraicLaw<Mesh>::
timeAdvance(const Vector& d_n) {
  _d_n_2.vec() = _d_n_1.vec();
  _d_n_1.vec() = _d_n.vec();
  _d_n.vec()   =  d_n;
}

//
//
//

//! Constructor
template <typename Mesh>  IndependentRing<Mesh>::
IndependentRing(const GetPot& dfile, const EntityFlag& marker):  
  DataSimplifiedStructure<Mesh>(dfile),
  DataTime(dfile,"solid/discretization"),
  _marker(marker),
  _dof(_mesh,_mesh.getRefFE()),
  _eta_nplus1(_dof.numTotalDof()),
  _w_nplus1(_dof.numTotalDof()),  
  _eta_nplus1old(_dof.numTotalDof()),
  _w_nplus1old(_dof.numTotalDof()), 
  _f(_dof.numTotalDof()), 
  _g(_dof.numTotalDof()),
  _d_n(_dof.numTotalDof() ),
  _d_n_1(_dof.numTotalDof() ),
  _d_n_2(_dof.numTotalDof() ) {
    
  a = _rho*_h;
  b = _E*_h/((1.0-_nu*_nu)*_R0*_R0);
  c = 2.0*a/_dt;
  d = 0.5*b + c/_dt;
 
}

//! Get initial condition
template <typename Mesh>  void IndependentRing<Mesh>::
initialize() {
  _eta_nplus1     = 0.0;
  _w_nplus1       = 0.0;  
  _eta_nplus1old  = 0.0;
  _w_nplus1old    = 0.0;
  _d_n_2.vec()    = 0.0;
  _d_n_1.vec()    = 0.0;
  _d_n.vec()      = 0.0;
}
 
//! Updating right hand side
template <typename Mesh>  
void IndependentRing<Mesh>::
timeAdvance(const Vector& d_n) { 
  _d_n_2.vec() = _d_n_1.vec();
  _d_n_1.vec() = _d_n.vec();
  _d_n.vec()   = d_n;

  // static right hand side displacement
  _f = (c/_dt+ -0.5*b)*_eta_nplus1 + c*_w_nplus1;
  // static right hand side velocity
  _g = -(2.0/_dt)*_eta_nplus1 -_w_nplus1;
}


//! Updates displacement
template <typename Mesh>  void IndependentRing<Mesh>::
iterate(const Vector& p) { 

  _eta_nplus1old = _eta_nplus1;
  _w_nplus1old   = _w_nplus1;

  // \eta^{n+1}
  _eta_nplus1 = (1.0/d)*(p + _f);
  // \w^{n+1}
  _w_nplus1 = (2.0/_dt)*_eta_nplus1 +_g;
  
  Real R,x,y; 

  UInt dim = _dof.numTotalDof();

  for (UInt i=0; i < dim; ++i) {

    if ( _mesh.point(i+1).marker() == _marker) {
      x = _mesh.point(i+1).x();    
      y = _mesh.point(i+1).y();    
 
      R = sqrt(x*x+y*y);
      
      _d_n.vec()[i]         = x * _eta_nplus1[i] / R;
      _d_n.vec()[i + dim]   = y * _eta_nplus1[i] / R;
    }
  }
}


//! Relaxation of the displacement
template <typename Mesh> void  IndependentRing<Mesh>::
relaxation(const Real& omega, Vector& d) {
  _eta_nplus1 = omega * _eta_nplus1 + (1.0-omega) * _eta_nplus1old;
  _w_nplus1   = omega * _w_nplus1   + (1.0-omega) * _w_nplus1old;

  Real R,x,y; 

  UInt dim = _dof.numTotalDof();
  
  for (UInt i=0; i < dim; ++i) {
    
    if ( _mesh.point(i+1).marker() == _marker) {
      x = _mesh.point(i+1).x();    
      y = _mesh.point(i+1).y();    
      
      R = sqrt(x*x+y*y);
      
      d[i]         = x * _eta_nplus1[i] / R;
      d[i + dim]   = y * _eta_nplus1[i] / R;
    }
  }

} 

//! Output
template <typename Mesh> void  IndependentRing<Mesh>::
showMe(ostream& c) const {
  DataSimplifiedStructure<Mesh>::showMe(c);
  DataTime::showMe(c);
}

#endif
