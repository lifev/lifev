/*!
  \file convDiffReactHandler.h
  \author M. Prosi
  \date 03/2004 
  \version 1.0

  \brief This file contains an abstract class for the Convection-Diffusion-Reaktion equation  solvers.  

*/

#ifndef _CONVDIFFREACTHANDLER_H_
#define _CONVDIFFREACTHANDLER_H_


#include "dataConvDiffReact.hpp"
#include "dataAztec.hpp"
#include "refFE.hpp"
#include "dof.hpp"
#include "lifeV.hpp"
#include "medit_wrtrs.hpp"
#include "bcCond.hpp"
#include "bdf.hpp"
#include "post_proc.hpp"
#include "openDX_wrtrs.hpp"
#include <cmath>
#include <sstream> 



/*! 
  \class convDiffReactHandler

  Abstract class which defines the general structure of a Convection-Diffusion-Reaction solver.
  For each new Convection-Diffusion-Reaction solver  we have to implement the corresponding 
  timeAdvance and an iterate methods 
  
*/

template <typename Mesh>
class ConvDiffReactHandler:
public DataConvDiffReact<Mesh> { 
 
 public:

  typedef Real (*Function)(const Real&, const Real&, const Real&, const Real&, const ID&);

  //! Constructor   
  /*!
    \param data_file GetPot data file
    \param refFE_c reference FE for the concentration
    \param Qr_c volumic quadrature rule for the concentration
    \param bdQr_c surface quadrature rule for the concentration
    \param BCh_c boundary conditions for the concentration
    \param ord_bdf order of the bdf time advancing scheme (default: Backward Euler) 
  */
  ConvDiffReactHandler(const GetPot& data_file,  const RefFE& refFE_c, 
		      const QuadRule& Qr_c, const QuadRule& bdQr_c, BC_Handler& BCh_c);
    
  //! Sets initial condition for the concentration (incremental approach): the initial time is t0, the time step dt
  void initialize(const Function& c0, Real t0, Real dt); 

 //! Sets initial condition for the concentration from file
  void initialize(const std::string & vname);

 //! Calculate the local elementsize of the elements for upwind 
  void calc_local_elemsize();

  //! Update the right  hand side  for time advancing   
  /*! 
    \param source volumic source  
    \param time present time
  */
  virtual void timeAdvance(const Function source, const Real& time) = 0; 

  //! Update convective term, bc treatment and solve the linearized cdr system
  virtual void iterate(const Real& time, PhysVectUnknown<Vector> & u) = 0;

  //! Returns the concentration vector
  ScalUnknown<Vector>& c();
  
  //! Returns the concentration Dof 
  const Dof& cDof() const;

  //! Returns the BDF Time Advancing stuff
  const Bdf& bdf() const;

  //! Returns the local element size
  ScalUnknown<Vector>& h();

  //! Do nothing destructor
  virtual ~ConvDiffReactHandler() {}

  
 protected:
  
  //! Reference FE for the concentration
  const RefFE& _refFE_c;

  //! The Dof object associated with the concentration
  Dof _dof_c;

  //! The number of total concentration dofs  
  UInt _dim_c;

  //! Quadrature rule for concentration volumic elementary computations
  const QuadRule& _Qr_c;
  
  //! Quadrature rule for concentration surface elementary computations
  const QuadRule& _bdQr_c;

  //! Current FE for the concentration c
  CurrentFE _fe_c;
  CurrentBdFE _feBd_c;

  //! The concenration
  ScalUnknown<Vector> _c;
   
  //! The BC handler
  BC_Handler& _BCh_c;

  // ! The BDF Time Advance Method 
  Bdf _bdf;

  //! The element size
  ScalUnknown<Vector> _h;

};



//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh> 
ConvDiffReactHandler<Mesh>::
ConvDiffReactHandler(const GetPot& data_file,  const RefFE& refFE_c, 
		    const QuadRule& Qr_c, const QuadRule& bdQr_c, BC_Handler& BCh_c):
     DataConvDiffReact<Mesh>(data_file),    
     _refFE_c(refFE_c),
     _dof_c(_mesh,_refFE_c),
     _dim_c(_dof_c.numTotalDof()),
     _Qr_c(Qr_c),
     _bdQr_c(bdQr_c),
     _fe_c(_refFE_c,_mesh.getGeoMap(),_Qr_c),
     _feBd_c(_refFE_c.boundaryFE(),_mesh.getGeoMap().boundaryMap(),_bdQr_c),
     _c(_dim_c),
     _h(_mesh.numVolumes()),
     _BCh_c(BCh_c),
     _bdf(_order_bdf) {}


// Returns the concentration
template<typename Mesh> ScalUnknown<Vector>& 
ConvDiffReactHandler<Mesh>::c() {
  return _c;
}


// Returns the local element size
template<typename Mesh> ScalUnknown<Vector>& 
ConvDiffReactHandler<Mesh>::h() {
  return _h;
}

// Returns the concentration Dof 
template<typename Mesh> const Dof& 
ConvDiffReactHandler<Mesh>::cDof() const {
  return _dof_c;
}


// Returns the BDF Time Advancing stuff
template<typename Mesh>  const Bdf& 
ConvDiffReactHandler<Mesh>::bdf() const {
  return _bdf;
}

// ! Initialize when  initial conditions concentration
template<typename Mesh> void 
ConvDiffReactHandler<Mesh>::initialize(const Function& c0, Real t0, Real dt) {
  
  _bdf.initialize_unk(c0,_mesh,_refFE_c,_fe_c,_dof_c,t0,dt,1);  
  _c=*(_bdf.unk().begin()); 

  _bdf.showMe();

}

// ! Initialize when initial values for the concentration are read from file
template<typename Mesh> void 
ConvDiffReactHandler<Mesh>::initialize(const std::string & vname) {
  

    std::fstream Resfile(vname.c_str(),ios::in | ios::binary);
    if (Resfile.fail()) {std::cerr<<" Error in initialize: File not found or locked"<<std::endl; abort();}
    Resfile.read((char*)&_c(1),_c.size()*sizeof(double));
    Resfile.close();

   _bdf.initialize_unk(_c);

   _bdf.showMe();

}


// ! Calculate the locale element size for upwind parameter
template<typename Mesh> void 
ConvDiffReactHandler<Mesh>::calc_local_elemsize() {

  UInt NVertex=4;
  Real x[4],y[4],z[4],dx,dy,dz;
  Real h=0.0,h_corr; 
  for (UInt ih=1;ih<=_mesh.numVolumes();ih++){
    for (UInt ih_c=0;ih_c<NVertex;ih_c++){
      x[ih_c]=_mesh.volumeList(ih).point(ih_c+1).x();
      y[ih_c]=_mesh.volumeList(ih).point(ih_c+1).y();
      z[ih_c]=_mesh.volumeList(ih).point(ih_c+1).z();
       }
    h=0;
    for (UInt ih_c=0;ih_c<NVertex;ih_c++){
     for (UInt jh_c=ih_c+1;jh_c<NVertex;jh_c++){
       dx= x[ih_c]-x[jh_c]; dy=y[ih_c]-y[jh_c]; dz=z[ih_c]-z[jh_c];
       h_corr=sqrt(dx*dx+dy*dy+dz*dz);
       if (h_corr>h) h=h_corr;
     }
    }
    _h(ih)=h;
  }
  

}

#endif
