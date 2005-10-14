/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Miguel A. Fernandez <miguel.fernandez@inria.fr>
      Date: 2005-10-09

 Copyright (C) 2005 INRIA

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
/**
  \file fluidToMaster.hpp
  \author M.A. Fernandez
  \date 10/2005

  \brief This file contains an interface class between a LifeV fluid solver
  and the masterFSI coupler of Jean-Fred's, used at INRIA and P6 
  to couple a fluid solver with an external (non-linear) solid solver
*/
#ifndef _LIFEVTOMASTER_H_
#define _LIFEVTOMASTER_H_

#include <life/lifecore/life.hpp>
#include <vector>
#include <life/lifefem/bcHandler.hpp>
#include <pvm3.h>

namespace LifeV
{

  template <typename Fluid>
  class FluidToMaster {

  public:
  
    FluidToMaster(Fluid& fluid, int masterID, int fsiStatus);

    void sdToMasterFSI();  

    void sdCoorToMasterFSI();

    void rvFromMasterFSI(int msg);

    int status() {return _fsiStatus;}
    
    void dep_fluid_interf(Vector& analytic_disp);
    void dep_fluid_interf_1(Vector& analytic_disp);
    void minus_rho_da_interf(Vector& mrhoda);
    void minus_rho_a_interf(Vector& mrhoa);  
  
  private:
    int _nbNodeFSI; 
 
    int _masterID;
    
    int _fsiStatus;

    int* _nodeFSI;   
    
    double* _force;
     
    double* _dep_master, *_dep1, *_dep2;
   
    Fluid& _fluid;
    
    Real _dt;  
  };
  
  template <typename Fluid>
  FluidToMaster<Fluid>::FluidToMaster(Fluid& fluid, int masterID, int fsiStatus):
    _nbNodeFSI(0), 
    _masterID(masterID),
    _fsiStatus(fsiStatus),
    _nodeFSI(0),   
    _force(0),
    _dep_master(0), 
    _dep1(0), 
    _dep2(0),
    _fluid(fluid), 
    _dt(0.0)
  { }

  template <typename Fluid>
  void FluidToMaster<Fluid>::sdToMasterFSI() {
    
    static int ipass=0;
    int nt = 3*_nbNodeFSI;
    double vref=1.;
    int dim_u = _fluid.uDof().numTotalDof();
    
    std::cout <<  "---> Fluid sends FORCE to coupler\n" <<endl;
    
    for (int k=0 ; k<_nbNodeFSI ; ++k)
      for(int l=0;l<3;l++) 
	_force[3*k+l]= _fluid.residual()[_nodeFSI[k]-1+l*dim_u]; 
      
    pvm_initsend(PvmDataDefault);
    pvm_pkint(&nt,1,1);
    pvm_pkdouble(_force,nt,1);
    if(!ipass){
      pvm_pkdouble(&vref,1,1);
      ipass = 1;
    }
    pvm_send(_masterID,120);
  }
  
  template <typename Fluid>
  void FluidToMaster<Fluid>::rvFromMasterFSI(int msg) {


    std::cout << "<-- Fluid receives from Master\n" << endl;
    switch(msg){
    case 140:
      {
	std::cout << " -- Fluid Receives Message 140\n" << endl;
	int nnifs,ntdlf;
	int* p_fnifs;

	pvm_recv(_masterID,140);
	pvm_upkint(&nnifs,1,1);
	p_fnifs = new int[nnifs];
	pvm_upkint(&ntdlf,1,1);
	pvm_upkint(p_fnifs,nnifs,1);
	
	_nbNodeFSI= nnifs;
	_nodeFSI = new int[nnifs];
	std::cout << "--> Nodes FSI:" << _nbNodeFSI << endl;

	for(int i=0;i<nnifs;i++) {
	  _nodeFSI[i] = p_fnifs[i]; 
	}
	

	_dep_master = new double[3*_nbNodeFSI]; // d^{n+1} 
	_dep1       = new double[3*_nbNodeFSI]; // d^{n}  
	_dep2       = new double[3*_nbNodeFSI]; // d^{n-1} 	
	_force      = new double[3*_nbNodeFSI]; // f^{n+1}

	for(int i=0; i<_nbNodeFSI; ++i) 
	  {
	    _dep_master[3*i ]=_dep_master[3*i+1]=_dep_master[3*i+2]=0.;
	    _dep1[3*i ]=_dep1[3*i+1]=_dep1[3*i+2]=0.;
	    _dep2[3*i ]=_dep2[3*i+1]=_dep2[3*i+2]=0.;
	    _force[3*i ]=_force[3*i+1]=_force[3*i+2]=0.;
	  }
	break;
      }
    case 110:
      { 
	std::cout << "Fluid Receives Message 110\n" << endl;
	int ntdlf;

	pvm_recv(_masterID,110);
	pvm_upkint(&_fsiStatus,1,1);
	std::cout << "--> " << _fsiStatus << endl << endl;

	
	if(_fsiStatus>=0){
	  
	  pvm_upkdouble(&_dt,1,1);
	  _fluid.setTimeStep(_dt);
	  pvm_upkint(&ntdlf,1,1);
	  
	  if(ntdlf != 3*_nbNodeFSI)
	    std::cout << "rv: ntdlf = " << ntdlf << " 3*nbNodeFSI =" 
		      << 3*_nbNodeFSI << std::endl;
	  
	  // if new time step
	  if (_fsiStatus == 1) 
	    for(int i=0; i<_nbNodeFSI; ++i) 
	      for (int l = 0; l < 3 ; ++l) 
		{
		  _dep2[3*i +l ]=_dep1[3*i +l ];
		  _dep1[3*i +l ]=_dep_master[3*i +l ];		  
		}
	  
	  pvm_upkdouble(_dep_master,3*_nbNodeFSI,1);
	  
	}
      }
    } 
  }
  
  
  
  template <typename Fluid>
  void FluidToMaster<Fluid>::dep_fluid_interf(Vector& analytic_disp)
  {
    int dim_mesh = _fluid.uDof().numTotalDof(); // warning: il faut metre le dof du maillage
    
    for(int i=0;i<_nbNodeFSI;i++)
      for(UInt icoord=0;icoord<3;icoord++)
	analytic_disp[_nodeFSI[i]- 1 + dim_mesh*icoord]=_dep_master[3*i + icoord];
    
  } 

  template <typename Fluid>
  void FluidToMaster<Fluid>::dep_fluid_interf_1(Vector& analytic_disp)
  {
    int dim_mesh = _fluid.uDof().numTotalDof(); // warning: il faut metre le dof du maillage
    
    for(int i=0;i<_nbNodeFSI;i++)
      for(UInt icoord=0;icoord<3;icoord++)
	analytic_disp[_nodeFSI[i]- 1 + dim_mesh*icoord]=_dep1[3*i + icoord];
    
  } 


  template <typename Fluid> 
  void FluidToMaster<Fluid>::minus_rho_da_interf(Vector& mrhoda)
  {
    int dim_fluid = _fluid.uDof().numTotalDof();
    double rhodti2 = _fluid.density()/(_dt*_dt);
    
    for(int i=0;i<_nbNodeFSI;i++)
      for(UInt icoord=0;icoord<3;icoord++)
	mrhoda[_nodeFSI[i]- 1 + dim_fluid*icoord]=  - rhodti2 * _dep_master[3*i + icoord];
  } 
  
  template <typename Fluid> 
  void FluidToMaster<Fluid>::minus_rho_a_interf(Vector& mrhoa)
  {
    int dim_fluid = _fluid.uDof().numTotalDof();
    double rhodti2 = _fluid.density()/(_dt*_dt);
    
    for(int i=0;i<_nbNodeFSI;i++)
      for(UInt icoord=0;icoord<3;icoord++)
	mrhoa[_nodeFSI[i]- 1 + dim_fluid*icoord]= - rhodti2 * ( _dep_master[3*i + icoord]  -  2*_dep1[3*i + icoord] + _dep2[3*i + icoord] );
  } 



  template <typename Fluid>
  void FluidToMaster<Fluid>::sdCoorToMasterFSI() 
  {
    
    double coor[_nbNodeFSI*3];
    int nt = 3*_nbNodeFSI;
    std::cout << "---> Fluid sends coor to coupler\n ";
    
    for(int i=0;i<_nbNodeFSI;i++)
      {
	coor[3*i  ] = _fluid.mesh().point(_nodeFSI[i]).x();
	coor[3*i+1] = _fluid.mesh().point(_nodeFSI[i]).y();
	coor[3*i+2] = _fluid.mesh().point(_nodeFSI[i]).z();
      }
    
    pvm_initsend(PvmDataDefault);
    pvm_pkint(&nt,1,1);
    pvm_pkdouble(coor,nt,1);
    pvm_send(_masterID,125);
  
  }
}
#endif //_LIFEVTOMASTER_H_
